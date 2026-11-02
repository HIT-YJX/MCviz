#include "kgraph.h"
#include <limits>
#include <stdexcept>

void KPartiteGraph::init(uint k) {
    clear();
    K = k;
    vertex_ids.resize(K);
    id_to_pos.resize(K);

    full_adj.resize(K);
}


uint KPartiteGraph::add_vertex(uint p, mcviz::VertexId id) {
    auto found = id_to_pos[p].find(id);
    if (found != id_to_pos[p].end()) return found->second;
    if (num_vertices_total == numeric_limits<int>::max()) throw overflow_error("Too many vertices");

    uint v = static_cast<uint>(vertex_ids[p].size());

    id_to_pos[p][id] = v;
    vertex_ids[p].push_back(id);
    full_adj[p].emplace_back(K);
    ++num_vertices_total;
    return v;
}



// Copy the shared graph into per-vertex, per-layer adjacency without changing IDs.
// Expected time O(nK+m*log(Delta+2)), space O(nK+m), with maximum degree Delta.
void KPartiteGraph::load(const mcviz::Graph& graph) {
    init(graph.K);
    for (uint p = 0; p < K; p++) {
        for (auto id : graph.vertex_ids[p]) add_vertex(p, id);
    }

    for (auto edge : graph.edges) {
        uint p = edge.u.layer, q = edge.v.layer;
        uint u = edge.u.pos, v = edge.v.pos;
        full_adj[p][u][q].push_back(v);

        full_adj[q][v][p].push_back(u);
    }
    num_edges = graph.edges.size();
    sort_all_adj();
}

void KPartiteGraph::sort_all_adj() {
    for (uint p = 0; p < K; p++) for (auto& vertex : full_adj[p]) {
        for (uint q = 0; q < K; q++) {
            auto& a = vertex[q];
            sort(a.begin(), a.end());
            a.erase(unique(a.begin(), a.end()), a.end());
        }
    }
}


bool KPartiteGraph::has_edge(uint p, uint u, uint q, uint v) const {
    if (p >= K || q >= K || p == q || u >= full_adj[p].size()) return false;
    const auto& a = full_adj[p][u][q];
    return binary_search(a.begin(), a.end(), v);
}

// Reset graph domains, adjacency structures, counters, and source identity.
void KPartiteGraph::clear() {
    K = 0;
    vertex_ids.clear();
    id_to_pos.clear();
    cross_adj.clear();
    nei.clear();
    full_adj.clear();

    num_vertices_total = 0;
    num_edges = 0;
    source_sha256.clear();
}

KPartiteGraph::~KPartiteGraph() = default;

