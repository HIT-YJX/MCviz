#include "representation.h"
#include "common/text.h"
#include "common/objective.h"
#include <algorithm>
#include <limits>

namespace mcviz {


// Attach graph vertex IDs to the hierarchy representation.
TreeRepresentation prepare_representation(const Graph& g, const Forest& f) {
    TreeRepresentation t;

    t.domain.offsets.push_back(0);
    for (const auto& ids : g.vertex_ids) {
        t.domain.ids.insert(t.domain.ids.end(), ids.begin(), ids.end());
        t.domain.offsets.push_back(t.domain.ids.size());
    }
    t.nodes = f.nodes; t.roots = f.roots;
    return t;
}

namespace {

class GapMinimum {
    std::size_t base = 1;
    std::vector<Score> a;
public:
    explicit GapMinimum(const std::vector<Score>& gaps) {
        while (base < gaps.size()) {
            base *= 2;
        }
        a.assign(2 * base, std::numeric_limits<Score>::max());
        std::copy(gaps.begin(), gaps.end(), a.begin() + base);

        for (auto i = base; i-- > 1;) a[i] = std::min(a[2 * i], a[2 * i + 1]);
    }



    Score get(std::size_t lo, std::size_t hi) const {
        Score v = std::numeric_limits<Score>::max();

        for (lo += base, hi += base; lo < hi; lo /= 2, hi /= 2) {
            if (lo & 1) v = std::min(v, a[lo++]);

            if (hi & 1) v = std::min(v, a[--hi]);
        }
        return v;
    }
};
}


// DFS makes subtrees contiguous; range minima recover each layer projection.
// Time O(n log(n+1)+K), output and auxiliary space O(n+K).
LinearRepresentation tree_to_linear(const TreeRepresentation& t) {
    LinearRepresentation l;

    l.domain = t.domain;
    const auto n = t.domain.ids.size(), K = t.domain.offsets.size() - 1;

    l.order.reserve(n); l.gaps.reserve(n); l.root_offsets.reserve(t.roots.size() + 1);
    std::vector<std::pair<std::size_t, Score>> stack;
    for (auto root : t.roots) {
        l.root_offsets.push_back(l.order.size());

        stack.emplace_back(root, 0);
        while (!stack.empty()) {
            auto [u, gap] = stack.back(); stack.pop_back();

            if (u < n) { l.order.push_back(u); l.gaps.push_back(gap); }
            else {
                const auto& node = t.nodes[u];
                for (std::size_t j = node.children.size(); j-- > 0;)
                    stack.emplace_back(node.children[j], j == 0 ? gap : node.height);
            }
        }
    }
    l.root_offsets.push_back(n);
    l.layer_order.resize(n); l.layer_gaps.assign(n, 0); l.layer_values.assign(n, 0);
    l.inverse_global.resize(n); l.inverse_layer.resize(n);
    std::vector<Layer> layers(n);
    std::vector<Score> values(n, 0);

    for (Layer i = 0; i < K; ++i)
        std::fill(layers.begin() + t.domain.offsets[i], layers.begin() + t.domain.offsets[i + 1], i);
    for (std::size_t u = n; u < t.nodes.size(); ++u)
        for (auto v : t.nodes[u].children) if (v < n) values[v] = t.nodes[u].height;
    auto cursor = t.domain.offsets;
    std::vector<std::size_t> previous(K, n);
    GapMinimum minimum(l.gaps);

    for (std::size_t i = 0; i < n; ++i) {
        auto v = l.order[i]; auto layer = layers[v]; auto at = cursor[layer]++;
        l.layer_order[at] = v; l.layer_values[at] = values[v];

        l.inverse_global[v] = i; l.inverse_layer[v] = at;
        if (previous[layer] != n) l.layer_gaps[at] = minimum.get(previous[layer] + 1, i + 1);
        previous[layer] = i;
    }
    return l;
}

// Decode shared gaps with a monotone stack, merging equal-height boundaries.
// Time and output/auxiliary space O(n+K); no original graph or hierarchy is needed.
TreeRepresentation linear_to_tree(const LinearRepresentation& l) {
    TreeRepresentation t;
    t.domain = l.domain;
    const auto n = l.order.size();
    t.nodes.resize(n); t.nodes.reserve(2 * n); t.roots.reserve(l.root_offsets.size() - 1);
    std::vector<HierarchyNode> stack;
    for (std::size_t r = 0; r + 1 < l.root_offsets.size(); ++r) {
        auto first = l.root_offsets[r], end = l.root_offsets[r + 1];
        std::size_t last = l.order[first];
        auto close = [&]() {
            auto node = std::move(stack.back()); stack.pop_back();
            node.children.push_back(last);

            last = t.nodes.size(); t.nodes.push_back(std::move(node));
        };
        for (std::size_t i = first + 1; i < end; ++i) {
            Score h = l.gaps[i];
            while (!stack.empty() && stack.back().height > h) close();
            if (!stack.empty() && stack.back().height == h) stack.back().children.push_back(last);
            else stack.push_back({h, {last}});
            last = l.order[i];
        }
        while (!stack.empty()) close();
        t.roots.push_back(last);
    }
    return t;
}

std::uint64_t representation_bytes(const TreeRepresentation& t) {
    Score b = checked_add(checked_mul(t.domain.offsets.size() + t.roots.size(), sizeof(std::size_t)),
                          checked_mul(t.domain.ids.size() + t.nodes.size(), sizeof(Score)));
    for (const auto& node : t.nodes) b = checked_add(b, checked_mul(node.children.size(), sizeof(std::size_t)));
    return b;
}

std::uint64_t representation_bytes(const LinearRepresentation& l) {
    Score b = checked_mul(l.domain.ids.size() + l.gaps.size() + l.layer_gaps.size() + l.layer_values.size(), sizeof(Score));
    return checked_add(b, checked_mul(l.domain.offsets.size() + l.order.size() + l.root_offsets.size() +
        l.layer_order.size() + l.inverse_global.size() + l.inverse_layer.size(), sizeof(std::size_t)));
}
}
