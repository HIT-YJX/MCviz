#ifndef BASELINE_KGRAPH_H
#define BASELINE_KGRAPH_H

#include "common.h"

class KPartiteGraph {
public:
    uint K;

    vector<vector<mcviz::VertexId>> vertex_ids;
    vector<unordered_map<mcviz::VertexId, uint>> id_to_pos;
    string source_sha256;

    vector<vector<vector<uint>>> cross_adj;

    vector<vector<vector<uint>>> nei;

    vector<vector<vector<vector<uint>>>> full_adj;

    uint num_vertices_total;
    ull num_edges;

    KPartiteGraph() : K(0), num_vertices_total(0), num_edges(0) {}

    void init(uint k);

    uint add_vertex(uint part, mcviz::VertexId id);
    void load(const mcviz::Graph& graph);

    void sort_all_adj();

    bool has_edge(uint up, uint ui, uint vp, uint vi) const;

    void clear();

    ~KPartiteGraph();

};

#endif
