#include "forest.h"
#include "common/text.h"
#include <algorithm>
#include <chrono>
#include <map>
#include <numeric>

namespace mcviz {
namespace {
struct DSU {
    std::vector<std::size_t> parent, size;

    explicit DSU(std::size_t n) : parent(n), size(n, 1) { std::iota(parent.begin(), parent.end(), 0); }

    std::size_t find(std::size_t u) {
        while (parent[u] != u) { parent[u] = parent[parent[u]]; u = parent[u]; }
        return u;
    }

    bool join(std::size_t u, std::size_t v) {
        u = find(u); v = find(v);
        if (u == v) return false;
        if (size[u] < size[v]) std::swap(u, v);
        parent[v] = u; size[u] += size[v];
        return true;
    }
};
}


std::size_t Forest::vertex_index(LocalVertex v) const {
    return offsets[v.layer] + v.pos;
}


// Kruskal preserves threshold components; equal weights merge as one batch.
// Cost is O(m log(m+1)+n+K) plus DSU/map work; final structure uses O(n+K) space.
Forest Forest::build(const Graph& g, const std::vector<WeightedEdge>& es, double* sort_ms) {
    Forest f;
    f.offsets.push_back(0);
    for (Layer i = 0; i < g.K; ++i) {
        for (VertexPos u = 0; u < g.vertex_ids[i].size(); ++u) f.vertices.push_back({i, u});
        f.offsets.push_back(f.vertices.size());
    }
    std::vector<std::size_t> order(es.size());
    std::iota(order.begin(), order.end(), 0);
    std::chrono::steady_clock::time_point sort_start;
    if (sort_ms) sort_start = std::chrono::steady_clock::now();
    std::stable_sort(order.begin(), order.end(), [&](auto a, auto b) { return es[a].weight > es[b].weight; });
    if (sort_ms) *sort_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - sort_start).count();
    DSU d(f.vertices.size());
    f.nodes.resize(f.vertices.size());
    std::vector<std::size_t> repr(f.vertices.size());
    std::iota(repr.begin(), repr.end(), 0);
    for (std::size_t start = 0; start < order.size();) {
        std::size_t end = start;
        Score h = es[order[start]].weight;
        std::vector<std::size_t> before;
        while (end < order.size() && es[order[end]].weight == h) {
            const auto& e = es[order[end]].edge;
            before.push_back(d.find(f.vertex_index(e.u)));
            before.push_back(d.find(f.vertex_index(e.v)));
            ++end;
        }
        std::sort(before.begin(), before.end());
        before.erase(std::unique(before.begin(), before.end()), before.end());
        for (std::size_t j = start; j < end; ++j) {
            auto id = order[j];
            if (d.join(f.vertex_index(es[id].edge.u), f.vertex_index(es[id].edge.v))) f.edge_ids.push_back(id);
        }
        std::map<std::size_t, std::vector<std::size_t>> groups;
        for (auto u : before) groups[d.find(u)].push_back(repr[u]);
        for (auto& entry : groups) if (entry.second.size() > 1) {
            repr[entry.first] = f.nodes.size();
            f.nodes.push_back({h, std::move(entry.second)});
        }
        start = end;
    }
    for (std::size_t u = 0; u < f.vertices.size(); ++u) if (d.find(u) == u) f.roots.push_back(repr[u]);
    return f;
}


Sets Forest::members(std::size_t node, unsigned K) const {
    Sets d(K);
    std::vector<std::size_t> stack{node};
    while (!stack.empty()) {
        auto u = stack.back(); stack.pop_back();
        if (u < vertices.size()) d[vertices[u].layer].push_back(vertices[u].pos);
        else for (auto v : nodes[u].children) stack.push_back(v);
    }
    for (auto& a : d) std::sort(a.begin(), a.end());
    return d;
}


std::vector<Sets> Forest::cut(Score threshold, unsigned K) const {
    std::vector<Sets> out;
    auto todo = roots;
    while (!todo.empty()) {
        auto u = todo.back(); todo.pop_back();
        if (u < vertices.size() || nodes[u].height >= threshold) out.push_back(members(u, K));
        else for (auto c : nodes[u].children) todo.push_back(c);
    }
    return out;
}
}
