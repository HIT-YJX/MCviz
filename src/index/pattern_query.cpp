#include "pattern_index.h"
#include "query_seed_core.h"

namespace mcviz {
// Cut the hierarchy, peel gamma cores, and obtain a bounded feasible seed before
// searching every surviving region. Search is exponential; region storage is O(n+H).
SearchResult PatternIndex::query(const Graph& g, const QueryRequest& request, RunInfo& run) const {
    const auto& spec = identity.spec;
    const auto q = effective_query(g, spec, request);
    const auto fixed = local_required(g, q);
    const auto& p = spec.pattern;
    run.algorithm = engine_name(engine);
    run.execution_path = query_execution_path();
    run.solver = "none";
    run.bound = bound;
    run.accelerated = p.connected();
    run.counters.clear();
    run.solver_calls = run.regions = 0;
    SearchResult best;
    best.sets.resize(g.K);
    SolverStats stats;
    auto search = options.search;
    CoreWorkspace workspace;
    search.core_backend = CoreBackend::reference;
    search.core_workspace = &workspace;
    search.reviewed_rules = true;
    search.backend = RegionBackend::pivot_ordered;
    bool complete = true;
    if (p.edges.empty()) {
        auto d = g.all_vertices();
        if (enough_capacity(d, q.gamma))
            best = {true, capacity_bound(p, spec.objective, d), std::move(d), Status::optimal, {}};
        run.solver = "all_vertices";
    } else if (!p.connected()) {
        auto r = solve_region(g, p, spec.objective, q.gamma, g.all_vertices(), fixed, {}, search);
        run.solver = r.backend;
        ++run.solver_calls;
        ++run.regions;
        stats = r.stats;
        best = std::move(r.best);
        complete = r.complete;
    } else {
        run.counters["hierarchy_region_query"] = 1;
        complete = seed_core_query(g, *this, q, fixed, search, best, stats, run);
    }
    accumulate_backend_counters(run.counters, stats);
    best.status = complete ? (best.found ? Status::optimal : Status::infeasible) : Status::incomplete;
    run.counters["search_nodes"] = stats.nodes;
    run.counters["search_leaves"] = stats.leaves;
    run.counters["core_removed"] = stats.core_removed;
    run.counters["bound_prunes"] = stats.bound_prunes;
    run.counters["PA01_calls"] = stats.pa01_calls;
    run.counters["PA04_calls"] = stats.pa04_calls;
    run.counters["PA09_calls"] = stats.pa09_calls;
    run.counters["PA09_prunes"] = stats.pa09_prunes;
    run.counters["PA08_tree_calls"] = stats.pa08_tree_calls;
    run.counters["PA08_core_calls"] = stats.pa08_core_calls;
    run.counters["PA08_prunes"] = stats.pa08_prunes;
    run.counters["PA08_budget_fallbacks"] = stats.pa08_budget_fallbacks;
    return best;
}
}
