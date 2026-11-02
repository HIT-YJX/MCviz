#pragma once

#include "query.h"
#include <utility>

namespace mcviz {
struct LayerPattern {
    unsigned K = 0;
    std::vector<std::pair<Layer, Layer>> edges;
    std::vector<std::vector<Layer>> neighbors;
    std::string selector, source;
    Layer anchor_layer = 0;
    std::vector<Layer> path_order_input;
    std::vector<Layer> path_order_internal;

    void normalize();
    bool adjacent(Layer i, Layer j) const;
    std::vector<std::vector<Layer>> components() const;
    bool connected() const;
    std::vector<Layer> star_centers() const;
    bool complete() const;
    std::string canonical() const;
    std::string hash() const;
};

std::vector<Layer> read_path_order(const std::string& path, unsigned K);
LayerPattern make_pattern(unsigned K, const std::string& selector,
                          Layer primary = 0, const std::string& order_file = {});
}
