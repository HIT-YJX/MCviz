#include "bcvizm_lp.h"
#include <algorithm>

namespace mcviz::baseline {
Score safe_edge_cap(const Graph& g, const LayerPattern& P, Objective obj,
                    LocalVertex u, LocalVertex v) {
    std::vector<Score> cap(g.K);
    for (Layer p = 0; p < g.K; p++) {
        cap[p] = g.vertex_ids[p].size();
        if (P.adjacent(u.layer, p)) {
            cap[p] = std::min<Score>(cap[p], g.neighbors(u.layer, u.pos, p).size());
        }
        if (P.adjacent(v.layer, p)) {
            cap[p] = std::min<Score>(cap[p], g.neighbors(v.layer, v.pos, p).size());
        }
    }
    return score(P, obj, cap);
}


// Use endpoint degree bounds to order vertices by maximum frontier weight.
// For fixed K, construction costs O(B+n^2+m) time and O(B+n+m_P) space;
// B is the graph identity encoding, and m_P counts required-pattern edges.
BcvizmIndex construct_bcvizm(const Graph& g, const BuildSpec& spec) {
    const auto& P = spec.pattern;
    BcvizmIndex index;
    index.identity.engine = "BCviz-M-LP";
    index.identity.value_kind = "safe-upper";
    index.identity.graph_hash = g.hash();
    index.identity.pattern_hash = P.hash();
    index.identity.spec = spec;

    std::vector<std::vector<unsigned>> ids(g.K);
    unsigned n = 0;
    for (Layer p = 0; p < g.K; p++) {
        ids[p].resize(g.vertex_ids[p].size());
        for (VertexPos v = 0; v < ids[p].size(); v++) {
            ids[p][v] = n++;
            index.domain.push_back({p, g.get_id(p, v)});
        }
    }
    std::vector<std::vector<unsigned>> adj(n);
    std::vector<std::vector<Score>> weights(n);
    for (auto e : g.edges) {
        if (!P.adjacent(e.u.layer, e.v.layer)) continue;
        unsigned u = ids[e.u.layer][e.u.pos], v = ids[e.v.layer][e.v.pos];
        Score w = safe_edge_cap(g, P, spec.objective, e.u, e.v);
        adj[u].push_back(v); weights[u].push_back(w);
        adj[v].push_back(u); weights[v].push_back(w);
    }

    std::vector<bool> done(n), touched(n);
    std::vector<Score> key(n);
    for (unsigned pos = 0; pos < n; pos++) {
        unsigned best = n;
        for (unsigned v = 0; v < n; v++) {
            if (done[v] || !touched[v]) continue;
            if (best == n || key[v] > key[best] ||
                (key[v] == key[best] && v < best)) best = v;
        }
        bool root = best == n;
        if (root) for (unsigned v = 0; v < n; v++) {
            if (!done[v] && (best == n || v < best)) best = v;
        }
        index.order.push_back({best, root ? 0 : key[best], root});
        done[best] = true;

        for (size_t z = 0; z < adj[best].size(); z++) {
            unsigned v = adj[best][z];
            if (!done[v]) {
                touched[v] = true;
                key[v] = std::max(key[v], weights[best][z]);
            }
        }
    }
    return index;
}

}
