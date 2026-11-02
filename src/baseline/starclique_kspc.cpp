#include "starclique.h"

// Enumerate active masks and vertex subsets as the exact general-pattern fallback.
// The public contract requires all K layers. For n vertices and maximum degree
// Delta, time is O(2^K*K + 2^n*(n^2*log(Delta+2)+K^2)); workspace is O(n+K).
StarClique StarCliqueFinder::find_max_pattern_clique_kspc_online(
        const vector<vector<bool>>& layer_pattern,
        uint min_layers) {
    best = StarClique(K);
    num_branches = 0;

    StarClique result(K);
    bool best_found = false;
    vector<bool> active(K, false);
    bool more_masks = true;
    while (more_masks) {
        vector<uint> layers;
        for (uint p = 0; p < K; p++) {
            if (active[p]) layers.push_back(p);
        }
        if (layers.size() >= min_layers) {
            vector<KpcVertex> vertices;
            for (uint p : layers) {
                for (uint v = 0; v < g.vertex_ids[p].size(); v++) {
                    vertices.push_back({p, v});
                }
            }
            vector<bool> selected(vertices.size(), false);
            bool more_sets = true;
            while (more_sets) {
                num_branches++;
                vector<vector<uint>> sets(K);
                for (size_t i = 0; i < vertices.size(); i++) {
                    if (selected[i]) sets[vertices[i].part].push_back(vertices[i].vex);
                }

                bool feasible = true;
                for (uint p : layers) {
                    if (sets[p].size() < max(1u, thresholds[p])) {
                        feasible = false;
                        break;
                    }
                }
                if (feasible) {
                    for (size_t i = 0; i < vertices.size() && feasible; i++) {
                        if (!selected[i]) continue;
                        const KpcVertex& a = vertices[i];
                        for (size_t j = i + 1; j < vertices.size(); j++) {
                            if (!selected[j]) continue;
                            const KpcVertex& b = vertices[j];
                            if (a.part != b.part && layer_pattern[a.part][b.part] &&
                                !g.has_edge(a.part, a.vex, b.part, b.vex)) {
                                feasible = false;
                                break;
                            }
                        }
                    }
                }
                if (feasible && contains_required(sets)) {
                    vector<uint> sizes(K, 0);
                    for (uint p : layers) sizes[p] = (uint)sets[p].size();
                    ull score = 0;
                    if (problem_type == MESC) {
                        for (uint i : layers) for (uint j : layers) {
                            if (i >= j || !layer_pattern[i][j]) continue;
                            ull term = 0;
                            if (sizes[i] != 0 && sizes[j] > numeric_limits<ull>::max() / sizes[i]) throw overflow_error("Pattern objective exceeds uint64");
                            term = (ull)sizes[i] * sizes[j];
                            if (term > numeric_limits<ull>::max() - score) throw overflow_error("Pattern objective exceeds uint64");
                            score += term;
                        }
                    } else if (problem_type == MVSC) {
                        for (uint p : layers) {
                            if (sizes[p] > numeric_limits<ull>::max() - score) throw overflow_error("Pattern objective exceeds uint64");
                            score += sizes[p];
                        }
                    } else {
                        score = UINT_MAX;
                        for (uint p : layers) score = min(score, (ull)sizes[p]);
                    }
                    if (!best_found || score > result.size) {
                        result.size = score;
                        result.sets.swap(sets);
                        best_found = true;
                    }
                }

                size_t bit = 0;
                while (bit < selected.size() && selected[bit]) selected[bit++] = false;
                if (bit == selected.size()) more_sets = false;
                else selected[bit] = true;
            }
        }

        size_t bit = 0;
        while (bit < active.size() && active[bit]) active[bit++] = false;
        if (bit == active.size()) more_masks = false;
        else active[bit] = true;
    }

    best = result;
    best.found = best_found;
    return best;
}

StarClique StarCliqueFinder::find_max_star_clique_kspc_online(uint min_layers) {
    return find_max_star_clique_kspc_merge(min_layers);
}
