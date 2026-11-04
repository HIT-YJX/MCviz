#include "star_build.h"
#include "legacy_star.h"
#include "common/text.h"
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <unordered_map>

namespace mcviz {
namespace {
// The anchor only resolves the two possible centers of a K=2 star.
Layer center_of(const LayerPattern& p) {
    auto centers = p.star_centers();
    require(!centers.empty(), "Star build requires a canonical spanning star");
    return std::find(centers.begin(), centers.end(), p.anchor_layer) == centers.end() ? centers[0]
                                                                                      : p.anchor_layer;
}


// Pack a typed secondary position; both public types are at most 32 bits.
std::uint64_t secondary_key(Layer layer, VertexPos pos) {
    return (std::uint64_t(layer) << 32) | pos;
}

// Budget accounting saturates instead of wrapping; it is not an objective value.
std::uint64_t saturated_sum(std::uint64_t a, std::uint64_t b) {
    return b > std::numeric_limits<std::uint64_t>::max() - a ? std::numeric_limits<std::uint64_t>::max()
                                                             : a + b;
}

// Intersect per-layer shared-neighbor candidates to form a necessary anchor domain.
// Time follows scanned wedges and intersections; workspace is O(n_primary).
class AnchorDomains {
    const Graph& g;
    Layer center;
    const Thresholds& beta;
    std::vector<unsigned> counts;
    std::vector<VertexPos> touched;

  public:
    AnchorDomains(const Graph& graph, Layer p, const Thresholds& lower)
        : g(graph), center(p), beta(lower), counts(g.vertex_ids[p].size()) {}

    bool domain(VertexPos u, std::vector<VertexPos>& result, std::uint64_t& work, std::uint64_t budget) {
        auto charge = [&](std::uint64_t amount) {
            work = saturated_sum(work, amount);
            return work <= budget;
        };
        Layer rare = center;
        std::uint64_t least = std::numeric_limits<std::uint64_t>::max();
        for (Layer j = 0; j < g.K; ++j)
            if (j != center) {
                const auto& neighbors = g.neighbors(center, u, j);
                if (!charge(neighbors.size())) return false;
                std::uint64_t cost = 0;
                for (auto v : neighbors)
                    cost = saturated_sum(cost, g.neighbors(j, v, center).size());
                if (cost < least) {
                    least = cost;
                    rare = j;
                }
            }
        if (rare == center || !charge(least)) return false;
        touched.clear();
        for (auto v : g.neighbors(center, u, rare))
            for (auto x : g.neighbors(rare, v, center)) {
                if (!counts[x]) touched.push_back(x);
                ++counts[x];
            }
        for (auto x : touched) {
            bool keep = counts[x] >= beta[rare];
            counts[x] = 0;
            for (Layer j = 0; j < g.K && keep; ++j)
                if (j != center && j != rare) {
                    const auto& a = g.neighbors(center, u, j);
                    const auto& b = g.neighbors(center, x, j);
                    if (!charge(saturated_sum(a.size(), b.size()))) {
                        for (auto v : touched)
                            counts[v] = 0;
                        return false;
                    }
                    std::size_t i = 0, k = 0;
                    Score common = 0;
                    while (i < a.size() && k < b.size() && common < beta[j]) {
                        if (a[i] < b[k])
                            ++i;
                        else if (b[k] < a[i])
                            ++k;
                        else {
                            ++i;
                            ++k;
                            ++common;
                        }
                    }
                    keep = common >= beta[j];
                }
            if (keep) result.push_back(x);
        }
        std::sort(result.begin(), result.end());
        return true;
    }
};

// Each anchor owns exactly its incident edges. Complete local enumeration gives
// their exact values; a partial trial is erased before independent exact fallback.
// Per-target improvement uses only that edge's witnessed lower bound, never global best.
void build_anchors(const Graph& active, const LayerPattern& p, Objective obj, Layer center,
                   const Thresholds& beta, const std::vector<unsigned>& thresholds,
                   legacy_star::ProblemType problem, std::uint64_t budget, StarWeightBuild& out,
                   std::uint64_t memory) {
    std::vector<std::vector<std::size_t>> incident(active.vertex_ids[center].size());
    std::vector<std::vector<std::size_t>> slots(active.K);
    for (Layer j = 0; j < active.K; ++j)
        if (j != center) slots[j].resize(active.vertex_ids[j].size());
    for (std::size_t id = 0; id < active.edges.size(); ++id) {
        auto e = active.edges[id];
        if (e.u.layer != center) std::swap(e.u, e.v);
        incident[e.u.pos].push_back(id);
    }
    AnchorDomains domains(active, center, beta);
    SolverOptions solver;
    solver.backend = RegionBackend::pivot_ordered;
    solver.pivot_memory_bytes = memory;
    for (VertexPos u = 0; u < incident.size(); ++u) {
        ++out.anchor_attempts;
        for (auto id : incident[u]) {
            auto e = active.edges[id];
            if (e.u.layer != center) std::swap(e.u, e.v);
            slots[e.v.layer][e.v.pos] = id;
        }
        std::vector<VertexPos> domain;
        std::uint64_t setup = incident[u].size(), internal = 0, callbacks = 0, publications = 0;
        bool refused = !budget || setup > budget || !domains.domain(u, domain, setup, budget);
        auto status = legacy_star::BudgetStatus::Partial;
        if (!refused && domain.size() < beta[center]) {
            ++out.anchor_empty_domains;
            out.anchor_work = saturated_sum(out.anchor_work, setup);
            continue; // This necessary domain proves every edge at this anchor infeasible.
        }
        if (!refused) {
            std::vector<std::size_t> columns;
            auto total_work = [&]() {
                return saturated_sum(setup, saturated_sum(internal, saturated_sum(callbacks, publications)));
            };
            auto improve = [&](auto upper, const auto& secondary, const auto& right) {
                if (refused) return true; // Let the budget poll trigger fallback; never prune on refusal.
                if (columns.empty()) {
                    callbacks = saturated_sum(callbacks, secondary.size());
                    columns.reserve(secondary.size());
                    for (auto [j, v] : secondary)
                        columns.push_back(slots[j][v]);
                }
                callbacks = saturated_sum(callbacks, right.size());
                if (total_work() > budget) {
                    refused = true;
                    return true;
                }
                for (auto v : right)
                    if (out.values[columns[v]] < upper) return true;
                return false;
            };
            auto publish = [&](const legacy_star::StarClique& witness) {
                if (refused) return;
                std::uint64_t credits = 0;
                for (Layer j = 0; j < active.K; ++j)
                    if (j != center) credits = saturated_sum(credits, witness.sets[j].size());
                callbacks = saturated_sum(callbacks, witness.sets[center].size());
                if (saturated_sum(total_work(), credits) > budget) {
                    refused = true;
                    return;
                }
                // u is adjacent to every local secondary column; every maximal witness includes u.
                if (std::find(witness.sets[center].begin(), witness.sets[center].end(), u) ==
                    witness.sets[center].end())
                    throw std::logic_error("Anchor witness omitted its universal primary vertex");
                publications += credits;
                ++out.anchor_witnesses;
                for (Layer j = 0; j < active.K; ++j)
                    if (j != center)
                        for (auto v : witness.sets[j])
                            out.values[slots[j][v]] = std::max<Score>(out.values[slots[j][v]], witness.size);
            };
            legacy_star::ExactCenterClosure enumeration(active, problem, thresholds, center, 0, 0, 0, false,
                                                        publish, &domain, u, improve);
            legacy_star::StarClique unused(active.K);
            status = enumeration.enumerate_budgeted(unused, [&](auto nodes, auto work) {
                (void)nodes;
                internal = work;
                return !refused && total_work() <= budget;
            });
            out.anchor_nodes = saturated_sum(out.anchor_nodes, enumeration.visited_nodes());
        }
        out.anchor_work = saturated_sum(
            out.anchor_work,
            saturated_sum(setup, saturated_sum(internal, saturated_sum(callbacks, publications))));
        out.anchor_publications = saturated_sum(out.anchor_publications, publications);
        if (status == legacy_star::BudgetStatus::Complete && !refused) {
            ++out.anchor_completed;
            continue;
        }
        out.anchor_budget_fallback += refused || status == legacy_star::BudgetStatus::Partial;
        out.anchor_memory_fallback += status == legacy_star::BudgetStatus::MemoryFallback;
        out.anchor_partial_publications = saturated_sum(out.anchor_partial_publications, publications);
        for (auto id : incident[u])
            out.values[id] = 0;
        for (auto id : incident[u]) {
            auto e = active.edges[id];
            if (e.u.layer != center) std::swap(e.u, e.v);
            auto r = solve_region(active, p, obj, beta, active.all_vertices(), {e.u, e.v}, {}, solver);
            require(r.complete, "Star anchor fallback search incomplete");
            out.values[id] = r.best.found ? r.best.score : 0;
            merge_solver_stats(out.solver_stats, r.stats);
            ++out.edge_calls;
        }
    }
}
} // namespace

// Every beta-feasible star survives the core. Materialize its required adjacencies once,
// retaining all core vertices (including any isolated vertices at beta-infeasible calls).
// Graph::from_records sorts the compact input: O((n+m_P) log(n+m_P)) time, O(n+m_P+K) space.
StarBuildView prepare_star_build(const Graph& g, const LayerPattern& p, const Sets& core) {
    (void)center_of(p);
    require(core.size() == g.K, "Star core layer mismatch");
    std::vector<std::vector<unsigned char>> member(g.K);
    std::vector<VertexRef> vertices;
    for (Layer i = 0; i < g.K; ++i) {
        member[i].resize(g.vertex_ids[i].size());
        for (auto u : core[i]) {
            require(u < member[i].size() && !member[i][u], "Invalid star core vertex");
            member[i][u] = 1;
            vertices.push_back({i, u});
        }
    }
    StarBuildView view;
    std::vector<std::pair<VertexRef, VertexRef>> edges;
    for (std::size_t id = 0; id < g.edges.size(); ++id) {
        const auto& e = g.edges[id];
        if (!p.adjacent(e.u.layer, e.v.layer) || !member[e.u.layer][e.u.pos] || !member[e.v.layer][e.v.pos])
            continue;
        edges.push_back({{e.u.layer, e.u.pos}, {e.v.layer, e.v.pos}});
        view.source_edge_ids.push_back(id);
    }
    view.graph = Graph::from_records(g.K, vertices, edges);
    require(view.graph.edges.size() == view.source_edge_ids.size(), "Star core edge count changed");
    // Compact positions preserve within-layer order; check the merge mapping explicitly.
    for (std::size_t i = 0; i < view.graph.edges.size(); ++i) {
        const auto& a = view.graph.edges[i];
        const auto& b = g.edges[view.source_edge_ids[i]];
        require(a.u.layer == b.u.layer && a.v.layer == b.v.layer &&
                    view.graph.get_id(a.u.layer, a.u.pos) == b.u.pos &&
                    view.graph.get_id(a.v.layer, a.v.pos) == b.v.pos,
                "Star core edge order changed");
    }
    return view;
}



// Credit each edge with its best maximal-star score; incomplete trials use exact fallback.
// Enumeration is exponential; workspace includes adjacency matrices and DFS frames.
StarWeightBuild exact_star_build(const Graph& active, const LayerPattern& p, Objective obj,
                                 const Thresholds& beta, unsigned primary_limit, std::uint64_t work_budget,
                                 bool anchors, std::uint64_t anchor_work_budget,
                                 std::uint64_t pivot_memory_bytes) {
    const Layer center = center_of(p);
    std::vector<unsigned> thresholds;
    for (auto x : beta) {
        require(x <= std::numeric_limits<unsigned>::max(), "Star threshold exceeds local capacity");
        thresholds.push_back(static_cast<unsigned>(x));
    }
    auto problem = static_cast<legacy_star::ProblemType>(static_cast<unsigned>(obj) + 1);
    StarWeightBuild out;
    out.values.resize(active.edges.size());
    if (anchors && primary_limit && active.vertex_ids[center].size() > primary_limit) {
        build_anchors(active, p, obj, center, beta, thresholds, problem, anchor_work_budget, out,
                      pivot_memory_bytes);
        return out;
    }
    if (primary_limit && active.vertex_ids[center].size() <= primary_limit) {
        out.shared_attempted = true;
        if (!work_budget)
            out.budget_fallback = true;
        else {
            std::vector<std::unordered_map<std::uint64_t, std::size_t>> targets(
                active.vertex_ids[center].size());
            for (std::size_t id = 0; id < active.edges.size(); ++id) {
                auto e = active.edges[id];
                if (e.u.layer != center) std::swap(e.u, e.v);
                require(e.u.layer == center, "Non-star edge in compact core");
                targets[e.u.pos].emplace(secondary_key(e.v.layer, e.v.pos), id);
            }
            // Budget before publication, so one large witness cannot overshoot the edge budget.
            bool publication_refused = false;
            auto publish = [&](const legacy_star::StarClique& witness) {
                if (publication_refused) return;
                std::uint64_t secondary = 0;
                for (Layer j = 0; j < active.K; ++j)
                    if (j != center) secondary = saturated_sum(secondary, witness.sets[j].size());
                const auto primary = std::uint64_t(witness.sets[center].size());
                const auto credits =
                    primary && secondary > std::numeric_limits<std::uint64_t>::max() / primary
                        ? std::numeric_limits<std::uint64_t>::max()
                        : primary * secondary;
                if (saturated_sum(out.shared_work, saturated_sum(out.shared_publications, credits)) >
                    work_budget) {
                    publication_refused = true;
                    return;
                }
                ++out.shared_witnesses;
                out.shared_publications += credits;
                for (auto u : witness.sets[center])
                    for (Layer j = 0; j < active.K; ++j)
                        if (j != center)
                            for (auto v : witness.sets[j]) {
                                auto at = targets[u].find(secondary_key(j, v));
                                require(at != targets[u].end(), "Shared star contains a missing edge");
                                out.values[at->second] =
                                    std::max<Score>(out.values[at->second], witness.size);
                            }
            };
            legacy_star::ExactCenterClosure enumeration(active, problem, thresholds, center, 0, 0, 0, true,
                                                        publish);
            legacy_star::StarClique unused(active.K);
            auto status = enumeration.enumerate_budgeted(unused, [&](auto nodes, auto work) {
                out.shared_nodes = nodes;
                out.shared_work = work;
                return !publication_refused && saturated_sum(work, out.shared_publications) <= work_budget;
            });
            out.shared_completed = status == legacy_star::BudgetStatus::Complete && !publication_refused;
            out.budget_fallback = status == legacy_star::BudgetStatus::Partial || publication_refused;
            out.memory_fallback = status == legacy_star::BudgetStatus::MemoryFallback;
            if (out.shared_completed) return out;
        }
        std::fill(out.values.begin(), out.values.end(), 0);
    }
    // Shared trials may retain their specialized enumeration; each independent
    // fallback edge uses the same production regional backend as queries.
    SolverOptions solver;
    solver.backend = RegionBackend::pivot_ordered;
    solver.pivot_memory_bytes = pivot_memory_bytes;
    for (std::size_t id = 0; id < active.edges.size(); ++id) {
        auto e = active.edges[id];
        if (e.u.layer != center) std::swap(e.u, e.v);
        auto r = solve_region(active, p, obj, beta, active.all_vertices(), {e.u, e.v}, {}, solver);
        require(r.complete, "Star edge fallback search incomplete");
        out.values[id] = r.best.found ? r.best.score : 0;
        merge_solver_stats(out.solver_stats, r.stats);
        ++out.edge_calls;
    }
    return out;
}
} // namespace mcviz
