#include "pattern_solver.h"
#include "common/text.h"
#include <algorithm>
#include <array>
#include <limits>

namespace mcviz {
std::vector<Layer> choose_cover(const LayerPattern& p, const Sets& d,
                               const Required& required, bool shape_rules) {
    std::vector<Score> w(p.K);
    Sets q(p.K);

    for (auto u : required) q.at(u.layer).push_back(u.pos);

    for (Layer i = 0; i < p.K; ++i) {
        std::sort(q[i].begin(), q[i].end());
        q[i].erase(std::unique(q[i].begin(), q[i].end()), q[i].end());

        w[i] = d[i].size() - (shape_rules ? q[i].size() : 0);
    }
    std::vector<bool> chosen(p.K, false);
    auto centers = p.star_centers();
    bool tree = shape_rules && p.edges.size() == p.K - 1 && p.connected();
    bool cycle = shape_rules && p.edges.size() == p.K && p.connected();

    for (const auto& ns : p.neighbors) if (ns.size() != 2) cycle = false;
    if (shape_rules && p.complete()) {
        auto leave = static_cast<Layer>(std::max_element(w.begin(), w.end()) - w.begin());
        std::fill(chosen.begin(), chosen.end(), true);

        chosen[leave] = false;
    } else if (shape_rules && !centers.empty()) {
        Layer center = centers[0];
        Score leaves = 0;
        for (Layer i = 0; i < p.K; ++i) if (i != center) leaves = checked_add(leaves, w[i]);
        if (w[center] <= leaves) chosen[center] = true;
        else for (Layer i = 0; i < p.K; ++i) chosen[i] = i != center;
    } else if (tree) {
        std::vector<Layer> parent(p.K, p.K), order{0};
        parent[0] = 0;
        for (std::size_t j = 0; j < order.size(); ++j) for (auto v : p.neighbors[order[j]]) {
            if (parent[v] == p.K) { parent[v] = order[j]; order.push_back(v); }
        }
        std::vector<std::array<Score, 2>> cost(p.K);
        for (std::size_t j = order.size(); j-- > 0;) {
            auto u = order[j]; cost[u] = {0, w[u]};
            for (auto v : p.neighbors[u]) if (parent[v] == u) {
                cost[u][0] = checked_add(cost[u][0], cost[v][1]);
                cost[u][1] = checked_add(cost[u][1], std::min(cost[v][0], cost[v][1]));
            }
        }
        for (auto u : order) chosen[u] = (u != 0 && !chosen[parent[u]]) || cost[u][1] <= cost[u][0];
    } else if (cycle) {
        std::vector<Layer> order{0};
        Layer prev = p.K, u = 0;
        while (order.size() < p.K) {
            auto next = p.neighbors[u][0] == prev ? p.neighbors[u][1] : p.neighbors[u][0];

            prev = u; u = next; order.push_back(u);
        }
        const Score inf = std::numeric_limits<Score>::max();

        Score best = inf;
        for (unsigned first = 0; first < 2; ++first) {
            std::array<Score, 2> cost{inf, inf};

            cost[first] = first ? w[order[0]] : 0;
            std::vector<std::array<unsigned, 2>> pred(p.K);
            for (Layer j = 1; j < p.K; ++j) {
                std::array<Score, 2> next{inf, inf};
                for (unsigned b = 0; b < 2; ++b) for (unsigned a = 0; a < 2; ++a) {
                    if ((!a && !b) || cost[a] == inf) continue;
                    Score value = checked_add(cost[a], b ? w[order[j]] : 0);
                    if (value < next[b]) { next[b] = value; pred[j][b] = a; }
                }
                cost = next;
            }
            for (unsigned last = 0; last < 2; ++last) if ((first || last) && cost[last] < best) {
                best = cost[last];
                unsigned b = last;
                for (std::size_t j = p.K; j-- > 0;) {
                    chosen[order[j]] = b != 0;

                    if (j) b = pred[j][b];
                }
            }
        }
    } else {
        for (auto e : p.edges) if (!chosen[e.first] && !chosen[e.second]) {
            chosen[w[e.first] <= w[e.second] ? e.first : e.second] = true;
        }
    }

    std::vector<Layer> x;
    for (Layer i = 0; i < p.K; ++i) if (chosen[i]) x.push_back(i);

    return x;
}
}
