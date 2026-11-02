#include "starclique.h"

vector<vector<bool>> StarCliqueFinder::kpc_layer_degree_prune(
        const vector<vector<bool>>& layer_pattern) const {
    vector<vector<bool>> keep(K);
    for (uint p = 0; p < K; p++) {
        keep[p].assign(g.vertex_ids[p].size(), true);
    }

    bool changed = true;
    while (changed) {
        changed = false;
        for (uint p = 0; p < K; p++) {
            for (uint v = 0; v < keep[p].size(); v++) {
                if (!keep[p][v]) continue;

                bool remove = false;
                for (uint q = 0; q < K && !remove; q++) {
                    if (q == p || !layer_pattern[p][q] || thresholds[q] == 0) continue;

                    uint kept_neighbors = 0;
                    for (uint nb : g.full_adj[p][v][q]) {
                        if (nb < keep[q].size() && keep[q][nb]) kept_neighbors++;
                    }
                    if (kept_neighbors < thresholds[q]) remove = true;
                }

                if (remove) {
                    keep[p][v] = false;
                    changed = true;
                }
            }
        }
    }

    return keep;
}

bool StarCliqueFinder::kpc_vertex_compatible_with_sets(
        const KpcVertex& v,
        const vector<vector<uint>>& current_sets,
        const vector<vector<bool>>& layer_pattern) const {
    for (uint q = 0; q < K; q++) {
        if (q == v.part || !layer_pattern[v.part][q]) continue;
        for (uint u : current_sets[q]) {
            if (!g.has_edge(v.part, v.vex, q, u)) return false;
        }
    }
    return true;
}

bool StarCliqueFinder::kpc_seed_compatible(
        const vector<vector<uint>>& seed_sets,
        uint part,
        uint vex,
        const vector<vector<bool>>& layer_pattern) const {
    KpcVertex v{part, vex};
    return kpc_vertex_compatible_with_sets(v, seed_sets, layer_pattern);
}

vector<StarCliqueFinder::KpcVertex> StarCliqueFinder::kpc_collect_candidates(
        const vector<vector<uint>>& current_sets,
        const vector<vector<bool>>& keep,
        const vector<uint>& seed_lower_bound,
        const vector<vector<bool>>& layer_pattern) const {
    vector<KpcVertex> candidates;

    for (uint p = 0; p < K; p++) {
        for (uint v = 0; v < keep[p].size(); v++) {
            if (!keep[p][v]) continue;
            if (v <= seed_lower_bound[p]) continue;

            bool already_selected = false;
            for (uint selected : current_sets[p]) {
                if (selected == v) {
                    already_selected = true;
                    break;
                }
            }
            if (already_selected) continue;

            KpcVertex candidate{p, v};
            if (kpc_vertex_compatible_with_sets(candidate, current_sets, layer_pattern)) {
                candidates.push_back(candidate);
            }
        }
    }

    sort(candidates.begin(), candidates.end(), [](const KpcVertex& a, const KpcVertex& b) {
        if (a.part != b.part) return a.part < b.part;
        return a.vex < b.vex;
    });
    return candidates;
}

// Expand each layer's minimum seed using larger compatible vertices only.
// Total time sums compatibility checks and candidate copies over visited states;
// the number of states can be exponential. Each recursion level stores O(n+K).
void StarCliqueFinder::kpc_expand(
        vector<vector<uint>>& current_sets,
        vector<KpcVertex>& candidates,
        const vector<vector<bool>>& keep,
        const vector<uint>& seed_lower_bound,
        const vector<vector<bool>>& layer_pattern,
        StarClique& best_result) {
    (void)keep;
    (void)seed_lower_bound;
    num_branches++;

    vector<uint> current_sizes(K, 0);
    vector<uint> candidate_sizes(K, 0);
    for (uint p = 0; p < K; p++) current_sizes[p] = current_sets[p].size();
    for (const auto& v : candidates) candidate_sizes[v.part]++;

    bool feasible_now = true;
    for (uint p = 0; p < K; p++) {
        if (current_sizes[p] < thresholds[p]) {
            feasible_now = false;
            break;
        }
    }
    if (feasible_now && contains_required(current_sets)) {
        ull cur_size = compute_pattern_size(current_sizes, layer_pattern);
        if (!best_result.found || cur_size > best_result.size) {
            best_result.size = cur_size;
            best_result.sets = current_sets;
            best_result.found = true;
        }
    }

    vector<uint> max_sizes(K, 0);
    for (uint p = 0; p < K; p++) {
        max_sizes[p] = current_sizes[p] + candidate_sizes[p];
        if (max_sizes[p] < thresholds[p]) return;
    }

    ull ub = compute_pattern_size(max_sizes, layer_pattern);
    if (best_result.found && ub <= best_result.size) return;

    for (size_t i = 0; i < candidates.size(); i++) {
        KpcVertex chosen = candidates[i];
        current_sets[chosen.part].push_back(chosen.vex);

        vector<KpcVertex> next_candidates;
        next_candidates.reserve(candidates.size() - i - 1);
        for (size_t j = i + 1; j < candidates.size(); j++) {
            const KpcVertex& cand = candidates[j];
            if (motif_compatible(chosen, cand, layer_pattern)) {
                next_candidates.push_back(cand);
            }
        }

        vector<uint> next_counts(K, 0);
        for (uint p = 0; p < K; p++) next_counts[p] = current_sets[p].size();
        for (const auto& cand : next_candidates) next_counts[cand.part]++;

        bool can_meet_thresholds = true;
        for (uint p = 0; p < K; p++) {
            if (next_counts[p] < thresholds[p]) {
                can_meet_thresholds = false;
                break;
            }
        }
        if (can_meet_thresholds) {
            kpc_expand(current_sets, next_candidates, keep, seed_lower_bound,
                       layer_pattern, best_result);
        }

        current_sets[chosen.part].pop_back();
        if (best_result.found && best_result.size >= ub) return;
    }
}

void StarCliqueFinder::kpc_enumerate_seeds(
        size_t depth,
        const vector<uint>& layer_order,
        const vector<vector<bool>>& keep,
        const vector<vector<bool>>& layer_pattern,
        vector<vector<uint>>& seed_sets,
        StarClique& best_result) {
    if (depth == layer_order.size()) {
        vector<uint> seed_lower_bound(K, 0);
        for (uint p = 0; p < K; p++) {
            if (seed_sets[p].empty()) return;
            seed_lower_bound[p] = seed_sets[p][0];
        }

        vector<KpcVertex> candidates =
            kpc_collect_candidates(seed_sets, keep, seed_lower_bound, layer_pattern);
        kpc_expand(seed_sets, candidates, keep, seed_lower_bound, layer_pattern, best_result);
        return;
    }

    uint p = layer_order[depth];
    for (uint v = 0; v < keep[p].size(); v++) {
        if (!keep[p][v]) continue;
        if (!kpc_seed_compatible(seed_sets, p, v, layer_pattern)) continue;

        seed_sets[p].push_back(v);
        kpc_enumerate_seeds(depth + 1, layer_order, keep, layer_pattern,
                            seed_sets, best_result);
        seed_sets[p].pop_back();
    }
}

StarClique StarCliqueFinder::find_max_pattern_clique_kpc_online(
        const vector<vector<bool>>& layer_pattern) {
    best = StarClique(K);
    num_branches = 0;


    vector<vector<bool>> keep = kpc_layer_degree_prune(layer_pattern);
    for (uint p = 0; p < K; p++) {
        bool has_vertex = false;
        for (bool kept : keep[p]) {
            if (kept) {
                has_vertex = true;
                break;
            }
        }
        if (!has_vertex) return best;
    }

    vector<uint> layer_order(K);
    for (uint p = 0; p < K; p++) layer_order[p] = p;
    sort(layer_order.begin(), layer_order.end(), [&](uint a, uint b) {
        uint ca = 0, cb = 0;
        for (bool kept : keep[a]) if (kept) ca++;
        for (bool kept : keep[b]) if (kept) cb++;
        if (ca != cb) return ca < cb;
        return a < b;
    });

    StarClique global_best(K);
    vector<vector<uint>> seed_sets(K);
    kpc_enumerate_seeds(0, layer_order, keep, layer_pattern, seed_sets, global_best);

    best = global_best;
    for (uint p = 0; p < K; p++) sort(best.sets[p].begin(), best.sets[p].end());
    return best;
}

