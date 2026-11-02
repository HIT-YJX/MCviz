#pragma once
#include "pattern_bounds.h"

namespace mcviz {
class CoreCapacity {
    const LayerPattern& p;
    Thresholds beta;
    Objective objective;
    std::vector<std::vector<std::size_t>> ids;
    std::vector<std::vector<Score>> deg;
    std::vector<Score> sizes, totals;
    std::optional<Layer> center;
    bool enabled = false;
public:
    CoreCapacity(const Graph& g, const BuildSpec& spec, const Sets& core);
    bool ready() const { return enabled; }
    Score bound(const GraphEdge& edge) const;
};
}
