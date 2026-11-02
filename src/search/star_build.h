#pragma once
#include "star_solver.h"
namespace mcviz {
struct StarBuildView { Graph graph; std::vector<std::size_t> source_edge_ids; };
StarBuildView prepare_star_build(const Graph& g, const LayerPattern& p, const Sets& core);
}
