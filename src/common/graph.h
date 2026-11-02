#pragma once

#include "query.h"
#include <map>

namespace mcviz {
struct GraphEdge { LocalVertex u, v; };


class Graph {
public:
    unsigned K = 0;
    std::vector<std::vector<VertexId>> vertex_ids;
    std::vector<GraphEdge> edges;
    std::string graph_source_sha256, vertices_source_sha256;

    static Graph load(unsigned K, const std::string& graph_file, const std::string& vertices_file = {});
    static Graph from_records(unsigned K, const std::vector<VertexRef>& vertices,
                              const std::vector<std::pair<VertexRef, VertexRef>>& edges);
    std::optional<VertexPos> find_pos(Layer i, VertexId id) const;
    VertexPos get_pos(Layer i, VertexId id) const;
    VertexId get_id(Layer i, VertexPos u) const;
    const std::vector<VertexPos>& neighbors(Layer i, VertexPos u, Layer j) const;
    bool has_edge(Layer i, VertexPos u, Layer j, VertexPos v) const;
    Sets all_vertices() const;
    std::string canonical() const;
    std::string hash() const;

private:
    std::vector<std::vector<std::map<Layer, std::vector<VertexPos>>>> adj_;
};
}
