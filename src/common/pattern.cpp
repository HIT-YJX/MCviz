#include "pattern.h"
#include "text.h"
#include <algorithm>
#include <numeric>
#include <sstream>

namespace mcviz {

void LayerPattern::normalize() {
    for (auto& e : edges) {
        if (e.first > e.second) std::swap(e.first, e.second);
    }

    std::sort(edges.begin(), edges.end());

    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    neighbors.assign(K, {});
    for (auto e : edges) {
        neighbors[e.first].push_back(e.second);
        neighbors[e.second].push_back(e.first);
    }
    for (auto& a : neighbors) std::sort(a.begin(), a.end());
}



// Check whether two layers share a pattern edge.
bool LayerPattern::adjacent(Layer i, Layer j) const {
    const auto& a = neighbors.at(i);

    return std::binary_search(a.begin(), a.end(), j);
}



// Find connected components of the layer pattern.
std::vector<std::vector<Layer>> LayerPattern::components() const {
    std::vector<std::vector<Layer>> out;
    std::vector<bool> seen(K, false);
    for (Layer i = 0; i < K; ++i) {
        if (seen[i]) continue;
        out.push_back({i});
        seen[i] = true;

        auto& c = out.back();

        for (std::size_t a = 0; a < c.size(); ++a) {
            for (Layer j : neighbors.at(c[a])) if (!seen[j]) {
                seen[j] = true;

                c.push_back(j);
            }
        }
    }
    return out;
}



bool LayerPattern::connected() const { return components().size() == 1; }



std::vector<Layer> LayerPattern::star_centers() const {
    std::vector<Layer> out;
    if (edges.size() != K - 1) return out;
    for (Layer i = 0; i < K; ++i) if (neighbors.at(i).size() == K - 1) out.push_back(i);

    return out;
}

// Check whether the pattern contains every layer pair.
bool LayerPattern::complete() const { return edges.size() == std::size_t(K) * (K - 1) / 2; }



// Serialize the layer pattern in canonical order.
std::string LayerPattern::canonical() const {
    std::ostringstream out;
    out << "layer-pattern-v1\n" << K << '\n';
    for (auto e : edges) out << e.first << ' ' << e.second << '\n';
    return out.str();
}


std::string LayerPattern::hash() const { return sha256(canonical()); }


std::vector<Layer> read_path_order(const std::string& path, unsigned K) {
    auto rows = data_lines(read_bytes(path));

    auto tokens = csv_fields(rows[0]);

    std::vector<Layer> order;
    for (const auto& t : tokens) {
        auto x = number(t, "1-based path-order layer", K);

        order.push_back(static_cast<Layer>(x - 1));
    }
    return order;
}

LayerPattern make_pattern(unsigned K, const std::string& selector, Layer primary,
                          const std::string& order_file) {
    LayerPattern p;
    p.K = K;
    p.selector = selector;
    p.anchor_layer = primary;
    if (selector == "complete") {
        for (Layer i = 0; i < K; ++i) for (Layer j = i + 1; j < K; ++j) p.edges.emplace_back(i, j);
    } else if (selector == "star") {
        for (Layer j = 0; j < K; ++j) if (j != primary) p.edges.emplace_back(primary, j);
    } else if (selector == "path" || selector == "cycle") {
        std::vector<Layer> order(K);
        std::iota(order.begin(), order.end(), 0);
        if (!order_file.empty()) order = read_path_order(order_file, K);
        if (selector == "path") {
            p.path_order_internal = order;
            for (Layer x : order) p.path_order_input.push_back(x + 1);
            p.source = order_file;
        }
        for (Layer j = 1; j < K; ++j) p.edges.emplace_back(order[j - 1], order[j]);
        if (selector == "cycle") p.edges.emplace_back(0, K - 1);
    } else if (selector != "empty") {
        p.source = selector.rfind("file:", 0) == 0 ? selector.substr(5) : selector;
        for (const auto& row : data_lines(read_bytes(p.source))) {
            auto t = fields(row);
            p.edges.emplace_back(static_cast<Layer>(number(t[0], "pattern layer", K - 1)),
                                 static_cast<Layer>(number(t[1], "pattern layer", K - 1)));
        }
    }
    p.normalize();
    return p;
}
}
