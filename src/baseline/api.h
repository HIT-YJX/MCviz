#pragma once

#include "../common/graph.h"
#include "../common/objective.h"

namespace mcviz::baseline {
SearchResult run(const Graph& graph, const LayerPattern& pattern, Objective objective,
                 const QueryRequest& request, const std::string& method,
                 unsigned min_layers = 0);
}
