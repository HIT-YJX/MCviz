#pragma once
#include "star_solver.h"

namespace mcviz {
// Temporary compact core graph. Its IDs are original local positions, never public IDs.
// source_edge_ids maps each compact graph edge back to the full original graph edge table.
struct StarBuildView {
    Graph graph;
    std::vector<std::size_t> source_edge_ids;
};
struct StarWeightBuild {
    SolverStats solver_stats;
    std::vector<Score> values;
    bool shared_attempted = false, shared_completed = false;
    bool budget_fallback = false, memory_fallback = false;
    std::uint64_t shared_nodes = 0, shared_witnesses = 0, shared_publications = 0;
    std::uint64_t shared_work = 0, edge_calls = 0;
    std::uint64_t anchor_attempts = 0, anchor_completed = 0, anchor_empty_domains = 0;
    std::uint64_t anchor_budget_fallback = 0, anchor_memory_fallback = 0;
    std::uint64_t anchor_nodes = 0, anchor_witnesses = 0, anchor_publications = 0;
    std::uint64_t anchor_work = 0, anchor_partial_publications = 0;
};

StarBuildView prepare_star_build(const Graph& g, const LayerPattern& p, const Sets& core);
StarWeightBuild exact_star_build(const Graph& active, const LayerPattern& p, Objective obj,
                                 const Thresholds& beta, unsigned primary_limit, std::uint64_t work_budget,
                                 bool anchors, std::uint64_t anchor_work_budget,
                                 std::uint64_t pivot_memory_bytes = 64ULL << 20);
} // namespace mcviz
