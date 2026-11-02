#include "pattern_index.h"
#include "build_filter.h"
#include "search/core_capacity.h"
#include "common/text.h"
#include "search/star_solver.h"
#include "search/star_build.h"
#include "search/star_support_build.h"
#include "search/pattern_rules.h"
#include "search/pattern_dynamic_build.h"
#include <algorithm>
#include <numeric>

namespace mcviz {

const char* engine_name(Engine e) { return e == Engine::exact ? "MCviz-LP" : "MCviz+-LP"; }


// Label the query route selected for this pattern.
std::string PatternIndex::query_execution_path() const {
    const auto& p = identity.spec.pattern;
    if (p.edges.empty()) return "empty_pattern";
    if (!p.connected()) return "generic_fallback";
    bool fast = !p.star_centers().empty() && options.solver == "auto";
    return fast ? "star_index" : "pattern_index_generic";
}

namespace {
double elapsed_ms(std::chrono::steady_clock::time_point start) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
}



void add_stats(SolverStats& a, const SolverStats& b) { merge_solver_stats(a,b); }

}


// Compute exact edge cohesion or admissible bounds, then build threshold components.
// Exact search may be exponential; storage includes O(n+m) state and solver workspace.
PatternIndex PatternIndex::construct(const Graph& g, const BuildSpec& spec, Engine engine,
                                     const IndexOptions& options) {
    if (engine == Engine::exact && spec.objective == Objective::edge &&
        !spec.pattern.edges.empty() && spec.pattern.connected()) {
        PatternIndex idx;
        idx.engine = engine;
        idx.options = options;
        idx.options.search.reviewed_rules = true;
        idx.options.search.backend = RegionBackend::pivot_ordered;
        idx.source_graph_sha256 = g.graph_source_sha256;
        idx.source_vertices_sha256 = g.vertices_source_sha256;

        idx.bound = "exact";

        idx.identity = {schema_version, contract_version, engine_name(engine), idx.bound,
                        g.hash(), spec.pattern.hash(), spec};
        idx.execution_path = idx.query_execution_path();
        idx.build_strategy = "general-dynamic-all";
        auto core = reduce(g, spec.pattern, spec.beta, g.all_vertices(), {}, CoreBackend::queue);
        idx.core_domains = std::move(core.domains);

        std::vector<Score> values(g.edges.size());

        if (core.possible) {
            auto result = build_pattern_dynamic(g, spec.pattern, spec.beta, idx.core_domains,
                options.build_threads, options.cohesion_cache_bytes, options.cohesion_work_budget,
                options.search.pivot_memory_bytes);
            values = std::move(result.values);
        }
        for (std::size_t id = 0; id < g.edges.size(); ++id)
            if (spec.pattern.adjacent(g.edges[id].u.layer, g.edges[id].v.layer))
                idx.weights.push_back({g.edges[id], values[id]});
        idx.forest = Forest::build(g, idx.weights);
        return idx;
    }
    bool star = !spec.pattern.star_centers().empty();

    bool fast = star && options.solver == "auto" && !options.search.deadline && !options.search.node_limit;
    PatternIndex idx;

    idx.engine = engine;

    idx.options = options;
    idx.options.search.reviewed_rules = options.solver == "auto";

    idx.options.search.backend = options.solver == "auto" ? RegionBackend::pivot_ordered : RegionBackend::cover;
    idx.source_graph_sha256 = g.graph_source_sha256;

    idx.source_vertices_sha256 = g.vertices_source_sha256;
    idx.bound = engine == Engine::exact ? "exact" : "auto";

    if (idx.bound == "auto") {
        idx.bound = fast ? "star-support" : "capacity";
        if (!star && options.solver == "auto" && spec.pattern.connected()) {
            if (spec.pattern.edges.size() == spec.pattern.K - 1) idx.bound = "tree-size";
            else if (spec.objective == Objective::balanced) idx.bound = "balanced-core";
        }
    }
    const bool support_build = engine == Engine::upper && fast && options.optimize_star_build && idx.bound == "star-support";
    idx.identity = {schema_version, contract_version, engine_name(engine), idx.bound, g.hash(), spec.pattern.hash(), spec};
    const auto& p = spec.pattern;
    idx.execution_path = idx.query_execution_path();
    const bool optimized = fast && options.optimize_star_build;
    idx.build_strategy = p.edges.empty() ? "empty-pattern" : (!p.connected() ? "generic-fallback" :
        (!fast ? "pattern-generic" : (engine == Engine::exact ? "reference-star-edge" : "reference-star-bound")));
    for (auto key : {"build_core_ms", "build_view_ms", "build_weights_ms", "build_forest_ms", "build_capacity_prepare_ms"})
        idx.build_timing[key] = 0;
    for (auto key : {"support_prepare_ms", "support_rows_ms", "support_weights_ms"})
        idx.build_timing[key] = 0;
    for (auto key : {"build_queue_core", "build_capacity_degrees", "build_filter_skipped_edge_calls", "build_prepared_adjacency", "build_core_vertices", "build_core_edges",
                     "build_shared_attempts", "build_shared_completed", "build_shared_budget_fallback",
                     "build_shared_memory_fallback", "build_shared_witnesses", "build_shared_nodes",
                     "build_shared_publications", "build_shared_work", "build_exact_edge_calls",
                     "build_anchor_attempts", "build_anchor_completed", "build_anchor_empty_domains",
                     "build_anchor_budget_fallback", "build_anchor_memory_fallback", "build_anchor_nodes",
                     "build_anchor_witnesses", "build_anchor_publications", "build_anchor_work",
                     "build_anchor_partial_publications"})
        idx.build_counters[key] = 0;
    for (auto key : {"support_requested_threads", "support_actual_threads", "support_reuse_requested",
                     "support_workers_allocated", "support_scratch_bytes_per_worker", "support_shared_bytes",
                     "support_rows_prepared", "support_edges_evaluated", "support_cache_peak_bytes",
                     "support_cache_budget_bytes", "support_fallback_rows", "support_thread_fallbacks",
                     "support_memory_fallbacks", "support_wedge_visits", "support_merge_pairs",
                     "support_pair_lookups", "support_scalar_edges"})
        idx.build_counters[key] = 0;
    if (support_build) {
        idx.build_counters["support_requested_threads"] = options.build_threads;
        idx.build_counters["support_reuse_requested"] = true;

        idx.build_counters["support_cache_budget_bytes"] = options.support_cache_bytes;
    }
    idx.build_counters["build_queue_core"] = options.build_core_backend == CoreBackend::queue;
    auto started = std::chrono::steady_clock::now();

    auto core = reduce(g, p, spec.beta, g.all_vertices(), {},
                       options.build_core_backend);
    idx.build_timing["build_core_ms"] = elapsed_ms(started);

    idx.core_domains = std::move(core.domains);
    idx.build_stats.core_removed = core.removed;
    idx.build_stats.pa01_calls = 1;
    for (const auto& d : idx.core_domains) idx.build_counters["build_core_vertices"] += d.size();
    if (!p.connected()) return idx;
    auto filtered = prepare_build_filter(idx, g, spec, core.possible);
    const auto& domains = filtered.domains;
    core.possible = filtered.possible;
    bool prepared = engine == Engine::upper && optimized && core.possible && idx.bound == "star-support";

    if (prepared) {
        started = std::chrono::steady_clock::now();

        auto view = prepare_star_build(g, p, domains);
        idx.build_timing["build_view_ms"] = elapsed_ms(started);
        idx.build_counters["build_prepared_adjacency"] = 1;
        idx.build_counters["build_core_edges"] = view.graph.edges.size();
        started = std::chrono::steady_clock::now();
        std::vector<Score> values;

            auto result = build_star_support(view.graph, p, spec.objective, spec.beta,
                                             options.build_threads, options.support_cache_bytes);
            idx.build_strategy = result.strategy;
            for (const auto& item : result.counters) idx.build_counters[item.first] = item.second;
            for (const auto& item : result.timing) idx.build_timing[item.first] = item.second;
            values = std::move(result.values);
        
        std::size_t active = 0;
        for (std::size_t id = 0; id < g.edges.size(); ++id) {
            const auto& e = g.edges[id];
            if (!p.adjacent(e.u.layer, e.v.layer)) continue;

            Score value = 0;
            if (active < view.source_edge_ids.size() && view.source_edge_ids[active] == id)
                value = values[active++];
            idx.weights.push_back({e, value});
        }
        idx.build_timing["build_weights_ms"] = elapsed_ms(started);
    } else {
        started = std::chrono::steady_clock::now();

        std::optional<CoreCapacity> capacities;
        if (core.possible && engine == Engine::upper && options.build_capacity_degrees &&
            (idx.bound == "capacity")) capacities.emplace(g, spec, domains);
        idx.build_counters["build_capacity_degrees"] = capacities && capacities->ready();

        idx.build_timing["build_capacity_prepare_ms"] = elapsed_ms(started);

        idx.build_timing["build_view_ms"] += idx.build_timing["build_capacity_prepare_ms"];

        started = std::chrono::steady_clock::now();
        for (std::size_t id = 0; id < g.edges.size(); ++id) {
            const auto& e = g.edges[id];
            if (!p.adjacent(e.u.layer, e.v.layer)) continue;

            if (!filtered.keep_edges.empty() && !filtered.keep_edges[id]) ++idx.build_counters["build_filter_skipped_edge_calls"];

            Required fixed{e.u, e.v};
            Score value = 0;
            if (core.possible && (filtered.keep_edges.empty() || filtered.keep_edges[id]) &&
                ((capacities && capacities->ready()) || contains_required(domains, fixed))) {
                if (engine == Engine::exact) {
                    auto r = solve_region(g, p, spec.objective, spec.beta, domains, fixed, {}, idx.options.search);
                    ++idx.build_solver_calls;
                    add_stats(idx.build_stats, r.stats);
                    if (r.best.found) value = r.best.score;
                } else if (idx.bound == "star-support") {
                    value = star_edge_bound(g, p, spec.objective, spec.beta, domains, e, idx.bound == "star-support");
                } else if (capacities && capacities->ready()) {
                    value = capacities->bound(e);
                } else {
                    auto d = condition_domains(g, p, domains, fixed);
                    if (d.possible && enough_capacity(d.domains, spec.beta)) {
                        value = capacity_bound(p, spec.objective, d.domains);
                        if (idx.bound == "tree-size" || idx.bound == "balanced-core") {
                            auto budget = idx.options.search.reviewed_work_budget;
                            auto ub = reviewed_upper(g, p, spec.objective, spec.beta, d.domains, fixed, budget);
                            idx.build_stats.pa08_tree_calls += ub.tree_used;
                            idx.build_stats.pa08_core_calls += ub.core_used;
                            idx.build_stats.pa08_budget_fallbacks += ub.budget_fallback;
                            value = ub.possible ? ub.upper : 0;
                        }
                    }
                }
            }
            idx.weights.push_back({e, value});
        }

        idx.build_timing["build_weights_ms"] = elapsed_ms(started);
        idx.build_counters["build_exact_edge_calls"] = idx.build_solver_calls;
    }
    apply_build_filter(idx, g, spec, filtered);

    started = std::chrono::steady_clock::now();
    idx.forest = Forest::build(g, idx.weights);
    idx.build_timing["build_forest_ms"] = elapsed_ms(started);
    return idx;
}


}
