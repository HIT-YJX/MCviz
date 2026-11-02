#pragma once
#include "common/graph.h"

namespace mcviz {
struct WeightedEdge { GraphEdge edge; Score weight = 0; };
struct HierarchyNode { Score height = 0; std::vector<std::size_t> children; };
struct Forest {
    std::vector<LocalVertex> vertices;
    std::vector<std::size_t> offsets;
    std::vector<std::size_t> edge_ids;
    std::vector<HierarchyNode> nodes;
    std::vector<std::size_t> roots;


    static Forest build(const Graph& g, const std::vector<WeightedEdge>& edges,
                        double* sort_ms = nullptr);
    Sets members(std::size_t node, unsigned K) const;
    std::vector<Sets> cut(Score threshold, unsigned K) const;
    std::size_t vertex_index(LocalVertex v) const;
};
}
