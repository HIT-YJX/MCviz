#include "starclique.h"
#include <functional>
#include <numeric>

// Enumerate maximal primary-secondary biclique arms, then intersect their primary sets.
// Rebuild leaf neighborhoods after each complete merge; the envelope also permits
// future primary shrinkage. Time sums arm-enumeration, envelope and merge work
// over exponential state spaces; workspace includes stored arms and recursion copies.
StarClique StarCliqueFinder::find_max_star_clique_kspc_merge(uint min_layers) {
    best = StarClique(K);
    num_branches = 0;

    const uint need_p = max(1u, thresholds[primary_part]);
    StarClique result(K);
    vector<uint> secondaries;

    struct Arm {
        vector<uint> primary;
    };

    auto get_arms = [&](uint q) {
        vector<Arm> arms;
        const uint np = (uint)g.vertex_ids[primary_part].size();
        const uint nq = (uint)g.vertex_ids[q].size();
        vector<uint> all_p(np), candidates(nq);
        iota(all_p.begin(), all_p.end(), 0);
        iota(candidates.begin(), candidates.end(), 0);
        auto common_count = [&](const vector<uint>& L, uint v) {
            const auto& neighbors = g.full_adj[q][v][primary_part];
            uint count = 0;
            if (L.size() <= neighbors.size()) {
                for (uint u : L) if (binary_search(neighbors.begin(), neighbors.end(), u)) count++;
            } else {
                for (uint u : neighbors) if (binary_search(L.begin(), L.end(), u)) count++;
            }
            return count;
        };
        function<void(const vector<uint>&, const vector<uint>&,
                      vector<uint>, vector<uint>)> enumerate;
        enumerate = [&](const vector<uint>& L, const vector<uint>& R,
                        vector<uint> P, vector<uint> Q) {
            vector<pair<uint, uint>> ranked;
            ranked.reserve(P.size());
            for (uint v : P) ranked.push_back({common_count(L, v), v});
            sort(ranked.begin(), ranked.end());
            P.clear();
            for (const auto& entry : ranked) P.push_back(entry.second);
            while (!P.empty()) {
                uint x = P.front();
                P.erase(P.begin());
                const auto& neighbors = g.full_adj[q][x][primary_part];
                vector<uint> next_l;
                if (L.size() <= neighbors.size()) {
                    for (uint u : L) {
                        if (binary_search(neighbors.begin(), neighbors.end(), u)) next_l.push_back(u);
                    }
                } else {
                    for (uint u : neighbors) {
                        if (binary_search(L.begin(), L.end(), u)) next_l.push_back(u);
                    }
                }
                if (next_l.size() < need_p) {
                    Q.push_back(x);
                    continue;
                }

                bool blocked = false;
                vector<uint> next_q;
                for (uint v : Q) {
                    uint degree = common_count(next_l, v);
                    if (degree == next_l.size()) { blocked = true; break; }
                    if (degree != 0) next_q.push_back(v);
                }
                if (blocked) { Q.push_back(x); continue; }

                vector<uint> next_r = R;
                next_r.push_back(x);
                vector<uint> next_p, dominated;
                for (uint v : P) {
                    uint degree = common_count(next_l, v);
                    if (degree == next_l.size()) {
                        next_r.push_back(v);
                        if (common_count(L, v) == next_l.size()) dominated.push_back(v);
                    } else if (degree != 0) {
                        next_p.push_back(v);
                    }
                }
                if (next_r.size() >= max(1u, thresholds[q])) {
                    arms.push_back({next_l});
                }
                if (!next_p.empty()) enumerate(next_l, next_r, move(next_p), move(next_q));

                for (uint v : dominated) {
                    auto it = find(P.begin(), P.end(), v);
                    if (it != P.end()) P.erase(it);
                    Q.push_back(v);
                }
                Q.push_back(x);
            }
        };
        enumerate(all_p, {}, move(candidates), {});
        return arms;
    };

    auto envelope_bound = [&](const vector<uint>& P, ull& ub) {
        if (P.size() < need_p) return false;
        vector<uint> primary = P;
        vector<vector<uint>> secondary(K);
        for (uint q : secondaries) {
            for (uint v = 0; v < g.vertex_ids[q].size(); v++) {
                uint deg = 0;
                for (uint u : P) if (g.has_edge(primary_part, u, q, v)) deg++;
                if (deg >= need_p) secondary[q].push_back(v);
            }
        }

        bool changed = true;
        while (changed) {
            changed = false;
            vector<bool> keep_p(primary.size(), true);
            for (size_t i = 0; i < primary.size(); i++) {
                uint u = primary[i];
                for (uint q : secondaries) {
                    uint deg = 0;
                    for (uint v : secondary[q]) {
                        if (g.has_edge(primary_part, u, q, v)) deg++;
                    }
                    if (deg < max(1u, thresholds[q])) {
                        keep_p[i] = false;
                        changed = true;
                        break;
                    }
                }
            }
            vector<uint> next_p;
            for (size_t i = 0; i < primary.size(); i++) {
                if (keep_p[i]) next_p.push_back(primary[i]);
            }
            if (next_p.size() != primary.size()) primary.swap(next_p);

            for (uint q : secondaries) {
                vector<uint> next_s;
                for (uint v : secondary[q]) {
                    uint deg = 0;
                    for (uint u : primary) {
                        if (g.has_edge(primary_part, u, q, v)) deg++;
                    }
                    if (deg >= need_p) next_s.push_back(v);
                    else changed = true;
                }
                if (next_s.size() != secondary[q].size()) secondary[q].swap(next_s);
            }
        }

        if (primary.size() < need_p) return false;
        vector<uint> sizes(K, 0);
        sizes[primary_part] = (uint)primary.size();
        for (uint q : secondaries) {
            sizes[q] = (uint)secondary[q].size();
            if (sizes[q] < max(1u, thresholds[q])) return false;
        }
        if (problem_type == MESC) {
            for (uint q : secondaries) {
                ull term = 0;
                if (sizes[primary_part] != 0 && sizes[q] > numeric_limits<ull>::max() / sizes[primary_part]) term = numeric_limits<ull>::max();
                else term = (ull)sizes[primary_part] * sizes[q];
                if (term > numeric_limits<ull>::max() - ub) ub = numeric_limits<ull>::max();
                else ub += term;
            }
        } else if (problem_type == MVSC) {
            for (uint q : secondaries) {
                if (sizes[q] > numeric_limits<ull>::max() - ub) ub = numeric_limits<ull>::max();
                else ub += sizes[q];
            }
            if (sizes[primary_part] > numeric_limits<ull>::max() - ub) ub = numeric_limits<ull>::max();
            else ub += sizes[primary_part];
        } else {
            ub = UINT_MAX;
            ub = min(ub, (ull)sizes[primary_part]);
            for (uint q : secondaries) ub = min(ub, (ull)sizes[q]);
        }
        return true;
    };

    auto evaluate = [&](const vector<uint>& P) {
        vector<vector<uint>> sets(K);
        sets[primary_part] = P;
        for (uint q : secondaries) {
            for (uint v = 0; v < g.vertex_ids[q].size(); v++) {
                bool common = true;
                for (uint u : P) {
                    if (!g.has_edge(primary_part, u, q, v)) {
                        common = false;
                        break;
                    }
                }
                if (common) sets[q].push_back(v);
            }
        }
        if (sets[primary_part].size() < need_p) return;
        for (uint q : secondaries) {
            if (sets[q].size() < max(1u, thresholds[q])) return;
            for (uint u : sets[primary_part]) {
                for (uint v : sets[q]) {
                    if (!g.has_edge(primary_part, u, q, v)) return;
                }
            }
        }

        vector<uint> sizes(K, 0);
        for (uint q : secondaries) sizes[q] = (uint)sets[q].size();
        sizes[primary_part] = (uint)sets[primary_part].size();
        ull score = 0;
        if (problem_type == MESC) {
            ull sum = 0;
            for (uint q : secondaries) {
                if (sizes[q] > numeric_limits<ull>::max() - sum) throw overflow_error("Star-clique objective sum exceeds uint64");
                sum += sizes[q];
            }
            if (sizes[primary_part] != 0 && sum > numeric_limits<ull>::max() / sizes[primary_part]) throw overflow_error("Star-clique MESC objective exceeds uint64");
            score = (ull)sizes[primary_part] * sum;
        } else if (problem_type == MVSC) {
            for (uint q : secondaries) {
                if (sizes[q] > numeric_limits<ull>::max() - score) throw overflow_error("Star-clique objective sum exceeds uint64");
                score += sizes[q];
            }
            if (sizes[primary_part] > numeric_limits<ull>::max() - score) throw overflow_error("Star-clique objective sum exceeds uint64");
            score += sizes[primary_part];
        } else {
            score = UINT_MAX;
            score = min(score, (ull)sizes[primary_part]);
            for (uint q : secondaries) score = min(score, (ull)sizes[q]);
        }
        if (!result.found || score > result.size) {
            result.size = score;
            result.sets.swap(sets);
            result.found = true;
        }
    };

    vector<vector<Arm>> arms_by_layer(K);
    function<void(size_t, const vector<uint>&)> merge;
    merge = [&](size_t depth, const vector<uint>& P) {
        num_branches++;
        ull ub = 0;
        if (!envelope_bound(P, ub)) return;
        if (result.found && ub <= result.size) return;
        if (depth == secondaries.size()) {
            evaluate(P);
            return;
        }

        uint q = secondaries[depth];
        set<vector<uint>> seen;
        for (const Arm& arm : arms_by_layer[q]) {
            vector<uint> next;
            set_intersection(P.begin(), P.end(), arm.primary.begin(), arm.primary.end(),
                             back_inserter(next));
            if (next.size() < need_p || !seen.insert(next).second) continue;
            merge(depth + 1, next);
        }
    };

    vector<uint> optional_layers;
    for (uint q = 0; q < K; q++) if (q != primary_part) optional_layers.push_back(q);
    if (min_layers == K) {
        secondaries = optional_layers;
        bool possible = true;
        for (uint q : secondaries) {
            arms_by_layer[q] = get_arms(q);
            if (arms_by_layer[q].empty()) { possible = false; break; }
        }
        if (possible) {
            vector<uint> all_p(g.vertex_ids[primary_part].size());
            iota(all_p.begin(), all_p.end(), 0);
            merge(0, all_p);
        }
        best = result;
        return best;
    }
    return best;
}
