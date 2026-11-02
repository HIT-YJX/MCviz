#pragma once
#include "pattern_solver.h"
namespace mcviz {
Score star_edge_bound(const Graph& g, const LayerPattern& p, Objective obj,
                      const Thresholds& beta, const Sets& core, GraphEdge edge, bool support);

class PreparedStarBounds {
public:
    PreparedStarBounds(const Graph& active, const LayerPattern& pattern, Objective objective,
                       const Thresholds& beta);
    Score bound(GraphEdge edge, bool support) const;

private:
    using Neighborhood = std::vector<VertexPos>;
    const LayerPattern& pattern_;
    const Thresholds& beta_;
    Objective objective_;
    Layer center_;
    unsigned layers_;
    std::size_t primary_count_;
    std::vector<const Neighborhood*> primary_adjacency_;
    std::vector<std::vector<const Neighborhood*>> secondary_adjacency_;
    const Neighborhood& adjacency(Layer layer, VertexPos vertex, Layer target) const;
};
}
