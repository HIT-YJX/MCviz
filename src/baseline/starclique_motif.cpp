#include "starclique.h"

ull StarCliqueFinder::compute_pattern_size(const vector<uint>& set_sizes,
                                            const vector<vector<bool>>& layer_pattern) const {
    if (problem_type == MESC) {
        ull total_edges = 0;
        for (uint i = 0; i < K; i++) {
            for (uint j = i + 1; j < K; j++) {
                if (layer_pattern[i][j]) {
                    ull term = (ull)set_sizes[i] * set_sizes[j];
                    if (term > numeric_limits<ull>::max() - total_edges) {
                        throw overflow_error("Pattern objective exceeds uint64");
                    }
                    total_edges += term;
                }
            }
        }
        return total_edges;
    }
    return compute_star_size(set_sizes, problem_type, primary_part);
}

bool StarCliqueFinder::motif_compatible(const MotifVertex& a,
                                         const MotifVertex& b,
                                         const vector<vector<bool>>& layer_pattern) const {
    if (a.part == b.part) return true;
    if (!layer_pattern[a.part][b.part]) return true;

    return g.has_edge(a.part, a.vex, b.part, b.vex);
}

// Pivot over candidate and excluded sets to enumerate maximal compatible sets.
// With n vertices, depth d and B visited states, a conservative time bound is
// O(B*(K^2+n^2*log(Delta+2)*log(n+1))); workspace is O(d*(n+K)).
// Delta is the maximum degree; same-layer and non-pattern pairs are compatible.
void StarCliqueFinder::motif_online_dfs(vector<vector<uint>>& current_sets,
                                         vector<MotifVertex>& candidates,
                                         vector<MotifVertex>& excluded,
                                         const vector<vector<bool>>& layer_pattern,
                                         StarClique& best_result) {
    num_branches++;

    vector<uint> current_sizes(K, 0);
    vector<uint> candidate_sizes(K, 0);
    for (uint p = 0; p < K; p++) current_sizes[p] = current_sets[p].size();
    for (const auto& v : candidates) candidate_sizes[v.part]++;

    for (uint p = 0; p < K; p++) {
        if (current_sizes[p] + candidate_sizes[p] < thresholds[p]) return;
    }

    if (candidates.empty()) {
        if (excluded.empty()) {
            bool feasible = true;
            for (uint p = 0; p < K; p++) {
                if (current_sizes[p] < thresholds[p]) {
                    feasible = false;
                    break;
                }
            }
            if (feasible && contains_required(current_sets)) {
                ull cur_size = compute_pattern_size(current_sizes, layer_pattern);
                if (!best_result.found || cur_size > best_result.size) {
                    best_result.size = cur_size;
                    best_result.sets = current_sets;
                    best_result.found = true;
                }
            }
        }
        return;
    }

    auto same_vertex = [](const MotifVertex& x, const MotifVertex& y) {
        return x.part == y.part && x.vex == y.vex;
    };

    bool has_pivot = false;
    MotifVertex pivot{0, 0};
    size_t best_cover = 0;
    auto consider_pivot = [&](const MotifVertex& p) {
        size_t cover = 0;
        for (const auto& v : candidates) {
            if (!same_vertex(p, v) && motif_compatible(p, v, layer_pattern)) cover++;
        }
        if (!has_pivot || cover > best_cover) {
            has_pivot = true;
            pivot = p;
            best_cover = cover;
        }
    };

    for (const auto& v : candidates) consider_pivot(v);
    for (const auto& v : excluded) consider_pivot(v);

    vector<MotifVertex> iteration_order;
    iteration_order.reserve(candidates.size());
    for (const auto& v : candidates) {
        if (!has_pivot || same_vertex(pivot, v) || !motif_compatible(pivot, v, layer_pattern)) {
            iteration_order.push_back(v);
        }
    }

    auto compatible_degree = [&](const MotifVertex& p) {
        size_t degree = 0;
        for (const auto& v : candidates) {
            if (!same_vertex(p, v) && motif_compatible(p, v, layer_pattern)) degree++;
        }
        return degree;
    };

    sort(iteration_order.begin(), iteration_order.end(),
         [&](const MotifVertex& a, const MotifVertex& b) {
             size_t da = compatible_degree(a);
             size_t db = compatible_degree(b);
             if (da != db) return da > db;
             if (a.part != b.part) return a.part < b.part;
             return a.vex < b.vex;
         });

    for (const auto& chosen : iteration_order) {
        auto pos = find_if(candidates.begin(), candidates.end(),
                           [&](const MotifVertex& v) {
                               return v.part == chosen.part && v.vex == chosen.vex;
                           });
        if (pos == candidates.end()) continue;

        candidates.erase(pos);
        current_sets[chosen.part].push_back(chosen.vex);

        vector<MotifVertex> next_candidates;
        next_candidates.reserve(candidates.size());
        for (const auto& v : candidates) {
            if (motif_compatible(chosen, v, layer_pattern)) next_candidates.push_back(v);
        }

        vector<MotifVertex> next_excluded;
        next_excluded.reserve(excluded.size());
        for (const auto& v : excluded) {
            if (motif_compatible(chosen, v, layer_pattern)) next_excluded.push_back(v);
        }

        motif_online_dfs(current_sets, next_candidates, next_excluded, layer_pattern, best_result);

        current_sets[chosen.part].pop_back();
        excluded.push_back(chosen);
    }
}

StarClique StarCliqueFinder::find_max_pattern_clique_motif_online(const vector<vector<bool>>& layer_pattern) {
    best = StarClique(K);
    num_branches = 0;


    vector<vector<bool>> keep(K);
    for (uint p = 0; p < K; p++) keep[p].assign(g.vertex_ids[p].size(), true);
    bool changed = true;
    while (changed) {
        changed = false;
        for (uint p = 0; p < K; p++) {
            for (uint u = 0; u < keep[p].size(); u++) {
                if (!keep[p][u]) continue;
                for (uint q = 0; q < K; q++) {
                    if (!layer_pattern[p][q]) continue;
                    bool supported = false;
                    for (uint v : g.full_adj[p][u][q]) {
                        if (keep[q][v]) { supported = true; break; }
                    }
                    if (!supported) {
                        keep[p][u] = false;
                        changed = true;
                        break;
                    }
                }
            }
        }
    }

    StarClique global_best(K);
    uint n_primary = g.vertex_ids[primary_part].size();

    for (uint anchor = 0; anchor < n_primary; anchor++) {
        if (!keep[primary_part][anchor]) continue;
        vector<vector<uint>> current_sets(K);
        current_sets[primary_part].push_back(anchor);

        vector<MotifVertex> candidates;
        vector<MotifVertex> excluded;

        for (uint u = 0; u < n_primary; u++) {
            if (!keep[primary_part][u]) continue;
            if (u < anchor) excluded.push_back({primary_part, u});
            else if (u > anchor) candidates.push_back({primary_part, u});
        }

        for (uint q = 0; q < K; q++) {
            if (q == primary_part) continue;
            if (layer_pattern[primary_part][q]) {
                for (uint v : g.full_adj[primary_part][anchor][q]) {
                    if (keep[q][v]) candidates.push_back({q, v});
                }
            } else {
                uint nq = g.vertex_ids[q].size();
                for (uint v = 0; v < nq; v++) {
                    if (keep[q][v]) candidates.push_back({q, v});
                }
            }
        }

        vector<uint> reachable(K, 0);
        for (uint p = 0; p < K; p++) reachable[p] = current_sets[p].size();
        for (const auto& v : candidates) reachable[v.part]++;

        bool enough = true;
        for (uint p = 0; p < K; p++) {
            if (reachable[p] < thresholds[p]) {
                enough = false;
                break;
            }
        }
        if (!enough) continue;

        motif_online_dfs(current_sets, candidates, excluded, layer_pattern, global_best);
    }

    best = global_best;
    for (uint p = 0; p < K; p++) sort(best.sets[p].begin(), best.sets[p].end());
    return best;
}

