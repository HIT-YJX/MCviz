#pragma once
#include "pattern_core.h"
#include <chrono>
#include <map>

namespace mcviz {
enum class RegionBackend { cover, pivot_ordered };
struct SolverStats {
    std::uint64_t nodes = 0, leaves = 0, core_removed = 0, bound_prunes = 0;
    std::uint64_t pa01_calls = 0, pa04_calls = 0, pa09_calls = 0, pa09_prunes = 0;
    std::uint64_t pa08_tree_calls = 0, pa08_core_calls = 0, pa08_prunes = 0, pa08_budget_fallbacks = 0;
    std::uint64_t cover_calls = 0, pivot_calls = 0, pivot_nodes = 0;
    std::uint64_t pivot_pivots = 0, pivot_branches = 0, pivot_prepare_bytes = 0, pivot_memory_fallbacks = 0;
};
struct SolverOptions {
    CoreBackend core_backend = CoreBackend::reference;
    CoreWorkspace* core_workspace = nullptr;
    std::optional<std::chrono::steady_clock::time_point> deadline;
    std::optional<std::uint64_t> node_limit;
    bool reviewed_rules = false;
    std::uint64_t reviewed_work_budget = 200000;

    bool research_skip_initial_core = false;

    RegionBackend backend = RegionBackend::cover;
    std::uint64_t pivot_memory_bytes = 64ULL << 20;
};
struct RegionResult {
    SearchResult best;
    bool complete = true;
    bool cutoff_applied = false;
    SolverStats stats;
    std::string backend = "cover_generic", fallback;
};
void merge_solver_stats(SolverStats& target, const SolverStats& value);
void accumulate_backend_counters(std::map<std::string,std::uint64_t>& counters, const SolverStats& stats);
std::vector<Layer> choose_cover(const LayerPattern& p, const Sets& domains,
                               const Required& required = {}, bool shape_rules = false);
RegionResult solve_region(const Graph& g, const LayerPattern& p, Objective obj,
                           const Thresholds& gamma, Sets domains, const Required& required = {},
                           std::optional<Score> cutoff = std::nullopt,
                           const SolverOptions& options = {});
}
