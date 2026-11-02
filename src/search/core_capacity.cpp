#include "core_capacity.h"
#include <algorithm>
#include <limits>

namespace mcviz {
CoreCapacity::CoreCapacity(const Graph& g, const BuildSpec& spec, const Sets& core)
    : p(spec.pattern), beta(spec.beta), objective(spec.objective) {
    for (auto [i, j] : p.edges) for (auto t : p.neighbors[i])
        if (p.adjacent(j, t)) return;
    const auto missing = std::numeric_limits<std::size_t>::max();
    ids.resize(g.K); sizes.resize(g.K);
    for (Layer i = 0; i < g.K; ++i) {
        ids[i].assign(g.vertex_ids[i].size(), missing);
        sizes[i] = core[i].size();
        for (auto u : core[i]) {
            ids[i][u] = deg.size();
            deg.emplace_back(p.neighbors[i].size(), 0);
        }
    }
    for (const auto& e : g.edges) {
        auto i = e.u.layer, j = e.v.layer;
        if (!p.adjacent(i, j)) continue;
        auto u = ids[i][e.u.pos], v = ids[j][e.v.pos];

        if (u == missing || v == missing) continue;
        auto a = std::lower_bound(p.neighbors[i].begin(), p.neighbors[i].end(), j) - p.neighbors[i].begin();

        auto b = std::lower_bound(p.neighbors[j].begin(), p.neighbors[j].end(), i) - p.neighbors[j].begin();
        ++deg[u][a]; ++deg[v][b];
    }

    auto centers = p.star_centers();
    if (!centers.empty() && objective == Objective::edge) {
        center = centers.front(); totals.resize(deg.size());
        for (auto u : core[*center]) {
            auto id = ids[*center][u];

            bool possible = true;
            for (std::size_t t = 0; t < deg[id].size(); ++t) {
                if (deg[id][t] < beta[p.neighbors[*center][t]]) possible = false;
                totals[id] = checked_add(totals[id], deg[id][t]);
            }

            if (!possible) totals[id] = 0;
        }
    }
    enabled = true;
}


Score CoreCapacity::bound(const GraphEdge& edge) const {
    const auto missing = std::numeric_limits<std::size_t>::max();

    auto u = edge.u, v = edge.v;
    if (ids[u.layer][u.pos] == missing || ids[v.layer][v.pos] == missing) return 0;
    if (center) {
        if (u.layer != *center) std::swap(u, v);
        const auto d = deg[ids[v.layer][v.pos]][0];
        if (d < beta[*center]) return 0;
        return checked_mul(d, totals[ids[u.layer][u.pos]]);
    }
    auto counts = sizes;
    for (auto q : {u, v}) {
        const auto& row = deg[ids[q.layer][q.pos]];
        for (std::size_t t = 0; t < row.size(); ++t) counts[p.neighbors[q.layer][t]] = row[t];
    }
    for (Layer i = 0; i < counts.size(); ++i) if (counts[i] < beta[i]) return 0;
    return score(p, objective, counts);
}
}
