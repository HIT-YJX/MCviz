#pragma once
#include "forest.h"

namespace mcviz {

struct RepresentationDomain {
    std::vector<std::size_t> offsets;
    std::vector<VertexId> ids;
};
struct TreeRepresentation {
    RepresentationDomain domain;
    std::vector<HierarchyNode> nodes;
    std::vector<std::size_t> roots;
};
struct LinearRepresentation {
    RepresentationDomain domain;
    std::vector<std::size_t> order, root_offsets;
    std::vector<Score> gaps;
    std::vector<std::size_t> layer_order, inverse_global, inverse_layer;
    std::vector<Score> layer_gaps, layer_values;
};
TreeRepresentation prepare_representation(const Graph& g, const Forest& f);
LinearRepresentation tree_to_linear(const TreeRepresentation& tree);
TreeRepresentation linear_to_tree(const LinearRepresentation& linear);
std::uint64_t representation_bytes(const TreeRepresentation& tree);
std::uint64_t representation_bytes(const LinearRepresentation& linear);
}
