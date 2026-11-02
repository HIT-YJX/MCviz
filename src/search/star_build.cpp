#include "star_build.h"
namespace mcviz {
StarBuildView prepare_star_build(const Graph& g, const LayerPattern& p, const Sets& core) {
    std::vector<std::vector<unsigned char>> member(g.K);
    std::vector<VertexRef> vertices;
    for (Layer i = 0; i < g.K; ++i) {
        member[i].resize(g.vertex_ids[i].size());
        for (auto u : core[i]) {
            member[i][u] = 1;
            vertices.push_back({i, u});
        }
    }
    StarBuildView view;
    std::vector<std::pair<VertexRef, VertexRef>> edges;
    for (std::size_t id = 0; id < g.edges.size(); ++id) {
        const auto& e = g.edges[id];
        if (!p.adjacent(e.u.layer, e.v.layer) || !member[e.u.layer][e.u.pos] ||
            !member[e.v.layer][e.v.pos]) continue;
        edges.push_back({{e.u.layer, e.u.pos}, {e.v.layer, e.v.pos}});
        view.source_edge_ids.push_back(id);
    }
    view.graph = Graph::from_records(g.K, vertices, edges);

    return view;
}
}
