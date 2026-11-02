#include "graph.h"
#include "text.h"
#include <algorithm>
#include <limits>
#include <numeric>
#include <sstream>
#include <tuple>

namespace mcviz {

// Normalize typed IDs and undirected edges before building sorted adjacency.
// Sorting costs O((n+m) log(n+m+1)); storage is O(n+m+K).
Graph Graph::from_records(unsigned K, const std::vector<VertexRef>& vertices,
                         const std::vector<std::pair<VertexRef, VertexRef>>& input) {
    Graph g;
    g.K = K;
    g.vertex_ids.resize(K);

    for (auto u : vertices) {
        g.vertex_ids[u.layer].push_back(u.id);
    }

    for (auto e : input) {
        g.vertex_ids[e.first.layer].push_back(e.first.id);
        g.vertex_ids[e.second.layer].push_back(e.second.id);
    }

    g.adj_.resize(K);
    for (Layer i = 0; i < K; ++i) {
        auto& ids = g.vertex_ids[i];
        std::sort(ids.begin(), ids.end());
        ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
        g.adj_[i].resize(ids.size());
    }
    std::vector<std::tuple<Layer, VertexPos, Layer, VertexPos>> es;
    for (auto e : input) {
        if (e.first.layer > e.second.layer) std::swap(e.first, e.second);
        es.emplace_back(e.first.layer, g.get_pos(e.first.layer, e.first.id),
                        e.second.layer, g.get_pos(e.second.layer, e.second.id));
    }
    std::sort(es.begin(), es.end());
    es.erase(std::unique(es.begin(), es.end()), es.end());
    for (auto e : es) {
        auto [i, u, j, v] = e;
        g.edges.push_back({{i, u}, {j, v}});
        g.adj_[i][u][j].push_back(v);
        g.adj_[j][v][i].push_back(u);
    }
    for (auto& layer : g.adj_) for (auto& u : layer) for (auto& a : u) {
        std::sort(a.second.begin(), a.second.end());
    }
    return g;
}


Graph Graph::load(unsigned K, const std::string& graph_file, const std::string& vertices_file) {
    const auto gb = read_bytes(graph_file);
    const auto vb = vertices_file.empty() ? std::string() : read_bytes(vertices_file);
    std::vector<VertexRef> vs;

    std::vector<std::pair<VertexRef, VertexRef>> es;
    for (const auto& row : data_lines(vb)) {
        auto t = fields(row);
        vs.push_back({static_cast<Layer>(number(t[0], "vertex layer", K - 1)), number(t[1], "vertex id")});
    }
    for (const auto& row : data_lines(gb)) {
        auto t = fields(row);

        es.push_back({{static_cast<Layer>(number(t[0], "edge layer", K - 1)), number(t[1], "vertex id")},
                      {static_cast<Layer>(number(t[2], "edge layer", K - 1)), number(t[3], "vertex id")}});
    }
    auto g = from_records(K, vs, es);
    g.graph_source_sha256 = sha256(gb);
    if (!vertices_file.empty()) g.vertices_source_sha256 = sha256(vb);
    return g;
}


std::optional<VertexPos> Graph::find_pos(Layer i, VertexId id) const {
    const auto& ids = vertex_ids[i];
    auto at = std::lower_bound(ids.begin(), ids.end(), id);
    if (at == ids.end() || *at != id) return std::nullopt;
    return static_cast<VertexPos>(at - ids.begin());
}


VertexPos Graph::get_pos(Layer i, VertexId id) const {
    auto pos = find_pos(i, id);
    return *pos;
}



VertexId Graph::get_id(Layer i, VertexPos u) const { return vertex_ids.at(i).at(u); }



const std::vector<VertexPos>& Graph::neighbors(Layer i, VertexPos u, Layer j) const {
    const auto& a = adj_.at(i).at(u);

    auto it = a.find(j);
    static const std::vector<VertexPos> empty;
    return it == a.end() ? empty : it->second;
}



bool Graph::has_edge(Layer i, VertexPos u, Layer j, VertexPos v) const {
    get_id(j, v);
    const auto& a = neighbors(i, u, j);
    return std::binary_search(a.begin(), a.end(), v);
}


Sets Graph::all_vertices() const {
    Sets out(K);
    for (Layer i = 0; i < K; ++i) {
        out[i].resize(vertex_ids[i].size());
        std::iota(out[i].begin(), out[i].end(), 0);
    }
    return out;
}



// Serialize the graph in a stable canonical order.
std::string Graph::canonical() const {
    std::ostringstream out;
    out << "layer-graph-v1\n" << K << '\n';
    for (Layer i = 0; i < K; ++i) for (auto id : vertex_ids[i]) out << "v " << i << ' ' << id << '\n';
    for (auto e : edges) out << "e " << e.u.layer << ' ' << get_id(e.u.layer, e.u.pos) << ' '
                            << e.v.layer << ' ' << get_id(e.v.layer, e.v.pos) << '\n';

    return out.str();
}



std::string Graph::hash() const { return sha256(canonical()); }
}
