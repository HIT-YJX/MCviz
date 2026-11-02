#include "pattern_solver.h"
#include "pattern_pivot.h"
#include "common/text.h"
#include "pattern_rules.h"
#include <algorithm>

namespace mcviz {
namespace {
struct State { Sets domains, selected; std::size_t layer = 0; };

bool exhausted(const SolverOptions& opt, std::uint64_t nodes) {
    return (opt.node_limit && nodes >= *opt.node_limit) ||
        (opt.deadline && std::chrono::steady_clock::now() >= *opt.deadline);
}
}

// Enumerate subsets in a vertex cover and saturate its independent complement.
// For N cover vertices, time is O(2^N * poly(n,K,|P|)); DFS copies use O((N+1)*(n+K)) space.
RegionResult solve_cover_region(const Graph& g, const LayerPattern& p, Objective obj,
                           const Thresholds& gamma, Sets domains, const Required& required,
                           std::optional<Score> cutoff, const SolverOptions& options) {
    RegionResult out;
    out.stats.cover_calls = 1;
    out.backend = options.reviewed_rules ? "cover_reviewed" : "cover_generic";
    out.cutoff_applied = cutoff.has_value();
    out.best.sets.resize(g.K);

    CoreResult core;
    if (options.research_skip_initial_core) {
        auto conditioned = condition_domains(g, p, std::move(domains), required);
        core.domains = std::move(conditioned.domains);
        core.possible = conditioned.possible && enough_capacity(core.domains, gamma);
    } else {
        core = reduce(g, p, gamma, std::move(domains), required, options.core_backend, options.core_workspace);

        out.stats.pa01_calls = 1;
    }

    out.stats.core_removed = core.removed;
    if (!core.possible) {
        out.best.status = cutoff ? Status::incomplete : Status::infeasible;
        return out;
    }
    auto cover = choose_cover(p, core.domains, required, options.reviewed_rules);
    out.stats.pa04_calls = options.reviewed_rules ? 1 : 0;
    std::uint64_t rule_budget = options.reviewed_work_budget;

    bool upper_eligible = !p.edges.empty() && (obj == Objective::balanced || (p.edges.size() == p.K - 1 && p.connected()));

    State root{std::move(core.domains), Sets(g.K), 0};
    for (auto q : required) root.selected[q.layer].push_back(q.pos);

    for (auto& c : root.selected) {
        std::sort(c.begin(), c.end());

        c.erase(std::unique(c.begin(), c.end()), c.end());
    }
    std::vector<State> stack;
    stack.push_back(std::move(root));

    while (!stack.empty()) {
        if (exhausted(options, out.stats.nodes)) { out.complete = false; break; }
        State st = std::move(stack.back());
        stack.pop_back();
        ++out.stats.nodes;
        if (!enough_capacity(st.domains, gamma)) continue;
        Score ub = capacity_bound(p, obj, st.domains);
        auto limit = out.best.found ? std::optional<Score>(out.best.score) : cutoff;
        if (limit && ub <= *limit) { ++out.stats.bound_prunes; continue; }

        if (options.reviewed_rules) {
            Required mandatory;
            for (Layer i = 0; i < g.K; ++i) for (auto u : st.selected[i]) mandatory.push_back({i, u});
            Thresholds lower = gamma;

            if (limit) {
                ++out.stats.pa09_calls;

                lower = improving_thresholds(p, obj, st.domains, gamma, st.selected, *limit);
                if (!enough_capacity(st.domains, lower)) { ++out.stats.pa09_prunes; continue; }
                if (lower != gamma) {
                    auto peeled = reduce(g, p, lower, std::move(st.domains), mandatory,
                                         options.core_backend, options.core_workspace);
                    ++out.stats.pa01_calls;

                    out.stats.core_removed += peeled.removed;

                    st.domains = std::move(peeled.domains);
                    if (!peeled.possible) { ++out.stats.pa09_prunes; continue; }
                    ub = capacity_bound(p, obj, st.domains);
                    if (ub <= *limit) { ++out.stats.pa09_prunes; continue; }
                }
            }
            if (upper_eligible) {
                if (rule_budget == 0) ++out.stats.pa08_budget_fallbacks;
                else {
                    auto upper = reviewed_upper(g, p, obj, lower, st.domains, mandatory, rule_budget,
                                                options.core_backend, options.core_workspace);
                    out.stats.pa08_tree_calls += upper.tree_used;
                    out.stats.pa08_core_calls += upper.core_used;
                    out.stats.pa08_budget_fallbacks += upper.budget_fallback;

                    if (!upper.possible || (limit && upper.upper <= *limit)) { ++out.stats.pa08_prunes; continue; }
                }
            }
        }

        std::optional<LocalVertex> next;
        while (st.layer < cover.size()) {
            Layer i = cover[st.layer];
            for (auto u : st.domains[i]) {
                if (!std::binary_search(st.selected[i].begin(), st.selected[i].end(), u)) {
                    next = LocalVertex{i, u};
                    break;
                }
            }

            if (next) break;
            ++st.layer;
        }
        if (!next) {
            ++out.stats.leaves;

            Score value = capacity_bound(p, obj, st.domains);

            SearchResult r{true, value, std::move(st.domains), Status::optimal, {}};

            out.best = std::move(r);
            continue;
        }

        auto u = *next;
        State skip = st;

        auto& skipped = skip.domains[u.layer];
        skipped.erase(std::lower_bound(skipped.begin(), skipped.end(), u.pos));
        stack.push_back(std::move(skip));
        auto& selected = st.selected[u.layer];
        selected.insert(std::lower_bound(selected.begin(), selected.end(), u.pos), u.pos);
        bool valid = true;

        for (auto j : p.neighbors[u.layer]) {
            auto& d = st.domains[j];
            d.erase(std::remove_if(d.begin(), d.end(), [&](VertexPos v) {
                return !g.has_edge(u.layer, u.pos, j, v);
            }), d.end());
            for (auto v : st.selected[j]) if (!std::binary_search(d.begin(), d.end(), v)) valid = false;
        }
        if (valid) stack.push_back(std::move(st));
    }
    if (!out.complete) out.best.status = Status::incomplete;
    else if (out.best.found) out.best.status = Status::optimal;
    else out.best.status = cutoff ? Status::incomplete : Status::infeasible;

    if (out.complete && cutoff && !out.best.found) out.best.message = "Region complete: no strict improvement above cutoff";
    return out;
}
}
