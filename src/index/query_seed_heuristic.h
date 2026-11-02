#pragma once
#include "pattern_index.h"

namespace mcviz {
// Deadline includes candidate preparation, greedy intersections and validation.
// Failure supplies no infeasibility certificate and changes no exact-search domain.
void seed_core_heuristic(const Graph&, const PatternIndex&, const QueryRequest&, const Required&,
    const Sets&, std::chrono::steady_clock::time_point, SearchResult&, RunInfo&);
} // namespace mcviz
