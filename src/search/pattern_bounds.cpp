#include "pattern_bounds.h"
#include "common/text.h"
#include <algorithm>

namespace mcviz {
bool contains_required(const Sets& d, const Required& required) {
    for (auto q : required) {
        const auto& c = d.at(q.layer);
        if (!std::binary_search(c.begin(), c.end(), q.pos)) return false;
    }
    return true;
}



bool enough_capacity(const Sets& d, const Thresholds& gamma) {
    for (Layer i = 0; i < d.size(); ++i) if (d[i].size() < gamma[i]) return false;
    return true;
}



DomainResult condition_domains(const Graph& g, const LayerPattern& p, Sets d,
                               const Required& required) {
    for (Layer i = 0; i < g.K; ++i) {
        auto& a = d[i];
        std::sort(a.begin(), a.end());
    }
    for (auto q : required) for (auto j : p.neighbors.at(q.layer)) {
        auto& a = d[j];
        a.erase(std::remove_if(a.begin(), a.end(), [&](VertexPos v) {
            return !g.has_edge(q.layer, q.pos, j, v);
        }), a.end());
    }
    bool possible = contains_required(d, required);
    return {std::move(d), possible};
}


Score capacity_bound(const LayerPattern& p, Objective obj, const Sets& d) {
    std::vector<Score> sizes;

    for (const auto& a : d) sizes.push_back(a.size());

    return score(p, obj, sizes);
}
}
