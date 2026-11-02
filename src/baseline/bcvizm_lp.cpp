#include "bcvizm_lp.h"
#include <algorithm>

namespace mcviz::baseline {
namespace {
struct Progress {
    BcvizmStats& stats;
    void state(bool star) {
        ++stats.stack_states;
        if (star) ++stats.star_states;
        else ++stats.generic_states;
    }
};

// Enumerate center subsets and saturate independent leaf layers with common neighbors.
// For c center and l leaf vertices, at most O(2^c) states perform filtering;
// auxiliary space is O((c+1)(K+l)), including common-neighbor copies.
void search_star_region(const Graph& g, const LayerPattern& P, Objective obj,
                        const Thresholds& tau, const Required& required,
                        Layer center, const std::vector<LocalVertex>& vertices,
                        const std::vector<unsigned>& order, SearchResult& best,
                        Progress& progress) {
    std::vector<VertexPos> primary;
    Sets common(g.K);
    for (unsigned v : order) {
        auto x = vertices[v];
        if (x.layer == center) primary.push_back(x.pos);
        else common[x.layer].push_back(x.pos);
    }
    auto& stats = progress.stats;
    ++stats.star_regions;
    stats.star_max_center_vertices = std::max<std::uint64_t>(stats.star_max_center_vertices, primary.size());

    stats.star_max_secondary_vertices = std::max<std::uint64_t>(
        stats.star_max_secondary_vertices, order.size() - primary.size());
    struct Frame { size_t at; unsigned char phase; };
    std::vector<Frame> stack{{0,0}};
    std::vector<Sets> common_stack;
    common_stack.push_back(std::move(common));
    std::vector<VertexPos> chosen;
    while (!stack.empty()) {
        Frame& frame = stack.back();
        if (frame.phase == 2) {
            chosen.pop_back();
            common_stack.pop_back();
            stack.pop_back();
            continue;
        }
        if (frame.phase == 1) {
            frame.phase = 2;
            VertexPos u = primary[frame.at];
            Sets next(g.K);
            for (Layer p = 0; p < g.K; ++p) if (p != center) {
                for (auto v : common_stack.back()[p]) {
                    if (g.has_edge(center, u, p, v)) next[p].push_back(v);
                }
            }
            chosen.push_back(u);

            common_stack.push_back(std::move(next));
            stack.push_back({frame.at + 1, 0});
            continue;
        }
        progress.state(true);

        const auto& current = common_stack.back();
        bool possible = chosen.size() + primary.size() - frame.at >= tau[center];
        for (Layer p = 0; p < g.K; ++p) if (p != center && current[p].size() < tau[p]) possible = false;
        for (auto q : required) if (q.layer != center) {
            const auto& candidates = current[q.layer];
            if (std::find(candidates.begin(), candidates.end(), q.pos) == candidates.end()) possible = false;
        }
        if (!possible) {
            ++stats.star_common_pruned;
            stack.pop_back();
            continue;
        }
        if (frame.at == primary.size()) {
            std::vector<Score> sizes(g.K);
            for (Layer p = 0; p < g.K; ++p) sizes[p] = p == center ? chosen.size() : current[p].size();
            ++stats.star_saturated_leaves;

            Score value = score(P, obj, sizes);
            if (!best.found || value > best.score) {
                best.found = true;

                best.score = value;
                best.sets = current;
                best.sets[center] = chosen;
            }

            stack.pop_back();
            continue;
        }

        ++stats.star_branch_states;

        frame.phase = 1;
        bool must_include = false;
        for (auto q : required) if (q.layer == center && q.pos == primary[frame.at]) must_include = true;
        if (!must_include) stack.push_back({frame.at + 1, 0});
    }
}

// Enumerate include/exclude states in saved order, pruning by layer capacities.
// For fixed K, r regional vertices and maximum degree d, time is
// O(2^r * r * log(d+2)); suffix counts and the stack use O(rK) space.
void search_region(const Graph& g, const LayerPattern& P, Objective obj,
                   const Thresholds& tau, const Required& required,
                   const std::vector<LocalVertex>& vertices,
                   const std::vector<unsigned>& order, SearchResult& best,
                   Progress& progress) {
    ++progress.stats.generic_regions;
    struct Frame { size_t at; unsigned char phase; };
    std::vector<Frame> stack{{0,0}};
    Sets chosen(g.K);
    std::vector<Score> sizes(g.K);

    std::vector<std::vector<Score>> suffix(order.size() + 1, std::vector<Score>(g.K));

    for (size_t at = order.size(); at > 0; at--) {
        suffix[at - 1] = suffix[at];
        suffix[at - 1][vertices[order[at - 1]].layer]++;
    }

    while (!stack.empty()) {
        Frame& frame = stack.back();
        if (frame.phase == 2) {
            Layer p = vertices[order[frame.at]].layer;
            chosen[p].pop_back();
            sizes[p]--;
            stack.pop_back();
            continue;
        }

        if (frame.phase == 1) {
            frame.phase = 2;
            const LocalVertex v = vertices[order[frame.at]];
            bool compatible = true;
            for (Layer p = 0; p < g.K && compatible; p++) {
                if (!P.adjacent(p, v.layer)) continue;

                for (auto u : chosen[p]) if (!g.has_edge(p, u, v.layer, v.pos)) {
                    compatible = false;
                    break;
                }
            }
            if (compatible) {
                chosen[v.layer].push_back(v.pos);
                sizes[v.layer]++;
                stack.push_back({frame.at + 1, 0});
            } else {
                frame.phase = 3;
            }

            continue;
        }
        if (frame.phase == 3) {
            stack.pop_back();
            continue;
        }

        progress.state(false);
        bool possible = true;
        std::vector<Score> cap(g.K);
        for (Layer p = 0; p < g.K; p++) {
            cap[p] = sizes[p] + suffix[frame.at][p];
            if (cap[p] < tau[p]) possible = false;
        }

        if (!possible || (best.found && score(P, obj, cap) <= best.score)) {
            stack.pop_back();

            continue;
        }

        if (frame.at == order.size()) {
            bool feasible = true;

            for (Layer p = 0; p < g.K; p++) if (sizes[p] < tau[p]) feasible = false;
            for (auto q : required) {
                const auto& C = chosen[q.layer];
                if (std::find(C.begin(), C.end(), q.pos) == C.end()) feasible = false;
            }
            if (feasible) {
                Score value = score(P, obj, sizes);
                if (!best.found || value > best.score) {
                    best.found = true;
                    best.score = value;
                    best.sets = chosen;
                }
            }
            stack.pop_back();
            continue;
        }
        ++progress.stats.generic_branch_states;

        frame.phase = 1;
        stack.push_back({frame.at + 1, 0});
    }
}
}

SearchResult query_bcvizm(const Graph& g, const BcvizmIndex& index,
                          const QueryRequest& request, BcvizmStats* stats) {
    BcvizmStats local_stats;
    if (!stats) stats = &local_stats;
    Progress progress{*stats};
    const auto& spec = index.identity.spec;
    const auto q = effective_query(g, spec, request);
    const auto required = local_required(g, q);
    const auto& P = spec.pattern;
    const auto centers = P.star_centers();
    Layer center = centers.empty() ? 0 : centers.front();
    if (std::find(centers.begin(), centers.end(), P.anchor_layer) != centers.end()) center = P.anchor_layer;
    SearchResult out;
    out.sets.resize(g.K);
    for (Layer p = 0; p < g.K; p++) if (q.gamma[p] > g.vertex_ids[p].size()) {
        out.status = Status::infeasible;

        return out;
    }

    std::vector<LocalVertex> vertices;
    for (Layer p = 0; p < g.K; p++) for (VertexPos v = 0; v < g.vertex_ids[p].size(); v++) {
        vertices.push_back({p,v});
    }
    if (P.edges.empty()) {
        out.found = true;

        out.sets = g.all_vertices();
        std::vector<Score> sizes(g.K);
        for (Layer p = 0; p < g.K; p++) sizes[p] = out.sets[p].size();
        out.score = score(P, spec.objective, sizes);
        out.status = Status::optimal;
        return out;
    }

    std::vector<std::vector<unsigned>> regions;

    if (!P.connected()) {
        regions.emplace_back();
        for (auto e : index.order) regions.back().push_back(e.vertex);
    } else {
        for (auto e : index.order) {
            if (regions.empty() || e.root || e.key == 0) {
                regions.emplace_back();
                if (stats) stats->zero_cuts++;
            }

            regions.back().push_back(e.vertex);
        }
    }

    for (const auto& region : regions) {
        if (stats) stats->regions++;
        std::vector<Score> cap(g.K);
        bool contains_required = required.empty();
        for (unsigned v : region) {
            cap[vertices[v].layer]++;
            if (!required.empty() && vertices[v].layer == required[0].layer &&
                vertices[v].pos == required[0].pos) contains_required = true;
        }
        bool possible = contains_required;
        for (Layer p = 0; p < g.K; p++) if (cap[p] < q.gamma[p]) possible = false;

        if (out.found && score(P, spec.objective, cap) <= out.score) possible = false;
        if (!possible) {
            if (stats) stats->capacity_skipped++;
            continue;
        }
        stats->active_region_vertices = region.size();
        if (!centers.empty()) {
            search_star_region(g, P, spec.objective, q.gamma, required, center, vertices, region, out, progress);
        } else {
            search_region(g, P, spec.objective, q.gamma, required, vertices, region, out, progress);
        }
    }
    out.status = out.found ? Status::optimal : Status::infeasible;

    return out;
}


SearchResult run_bcvizm(const Graph& g, const LayerPattern& P, Objective obj,
                        const QueryRequest& request, BcvizmStats* stats) {
    BuildSpec spec{P, obj, request.gamma};
    auto index = construct_bcvizm(g, spec);

    return query_bcvizm(g, index, request, stats);
}
}
