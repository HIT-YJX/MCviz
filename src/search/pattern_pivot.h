#pragma once
#include "pattern_solver.h"

namespace mcviz {
RegionResult solve_pivot_region(const Graph& g, const LayerPattern& p, Objective obj,
    const Thresholds& gamma, Sets domains, const Required& fixed,
    std::optional<Score> cutoff, const SolverOptions& options);
RegionResult solve_cover_region(const Graph& g, const LayerPattern& p, Objective obj,
    const Thresholds& gamma, Sets domains, const Required& fixed,
    std::optional<Score> cutoff, const SolverOptions& options);
}
