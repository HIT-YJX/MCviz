#ifndef BASELINE_STARCLIQUE_H
#define BASELINE_STARCLIQUE_H

#include "common.h"
#include "kgraph.h"
#include <limits>
#include <stdexcept>
#include <climits>
#include <optional>

class StarCliqueFinder {
public:
    const KPartiteGraph& g;
    uint K;
    ProblemType problem_type;
    uint primary_part;
    vector<vector<bool>> layer_pattern;
    optional<mcviz::LocalVertex> required;
    StarClique best;

    vector<uint> thresholds;

    ull num_branches;

    StarCliqueFinder(const KPartiteGraph& graph, ProblemType pt, const vector<uint>& thresh, uint pp = 0)
        : g(graph), K(graph.K), problem_type(pt), primary_part(pp), thresholds(thresh), num_branches(0) {
        best = StarClique(K);
        layer_pattern.assign(K, vector<bool>(K, false));
        for (uint q = 0; q < K; q++) {
            if (q == primary_part) continue;
            layer_pattern[primary_part][q] = true;
            layer_pattern[q][primary_part] = true;
        }
    }

    bool contains_required(const vector<vector<uint>>& sets) const;

    StarClique find_max_pattern_clique_motif_online(const vector<vector<bool>>& layer_pattern);

    StarClique find_max_pattern_clique_kpc_online(const vector<vector<bool>>& layer_pattern);

    StarClique find_max_pattern_clique_kspc_online(const vector<vector<bool>>& layer_pattern,
                                                   uint min_layers);
    StarClique find_max_star_clique_kspc_online(uint min_layers);

    ull compute_pattern_size(const vector<uint>& set_sizes,
                             const vector<vector<bool>>& layer_pattern) const;

private:
    struct MotifVertex {
        uint part;
        uint vex;
    };

    using KpcVertex = MotifVertex;

    bool motif_compatible(const MotifVertex& a,
                          const MotifVertex& b,
                          const vector<vector<bool>>& layer_pattern) const;
    void motif_online_dfs(vector<vector<uint>>& current_sets,
                          vector<MotifVertex>& candidates,
                          vector<MotifVertex>& excluded,
                          const vector<vector<bool>>& layer_pattern,
                          StarClique& best_result);

    vector<vector<bool>> kpc_layer_degree_prune(const vector<vector<bool>>& layer_pattern) const;
    bool kpc_vertex_compatible_with_sets(const KpcVertex& v,
                                         const vector<vector<uint>>& current_sets,
                                         const vector<vector<bool>>& layer_pattern) const;
    bool kpc_seed_compatible(const vector<vector<uint>>& seed_sets,
                             uint part,
                             uint vex,
                             const vector<vector<bool>>& layer_pattern) const;
    vector<KpcVertex> kpc_collect_candidates(const vector<vector<uint>>& current_sets,
                                             const vector<vector<bool>>& keep,
                                             const vector<uint>& seed_lower_bound,
                                             const vector<vector<bool>>& layer_pattern) const;
    void kpc_enumerate_seeds(size_t depth,
                             const vector<uint>& layer_order,
                             const vector<vector<bool>>& keep,
                             const vector<vector<bool>>& layer_pattern,
                             vector<vector<uint>>& seed_sets,
                             StarClique& best_result);
    void kpc_expand(vector<vector<uint>>& current_sets,
                    vector<KpcVertex>& candidates,
                    const vector<vector<bool>>& keep,
                    const vector<uint>& seed_lower_bound,
                    const vector<vector<bool>>& layer_pattern,
                    StarClique& best_result);

    StarClique find_max_star_clique_kspc_merge(uint min_layers);

};

inline ull compute_star_size(const vector<uint>& sizes, ProblemType pt, uint primary_part = 0) {
    if (sizes.empty()) return 0;
    if (pt == MBSC) return *min_element(sizes.begin(), sizes.end());

    const ull limit = numeric_limits<ull>::max();
    ull sum = 0;
    for (size_t i = 0; i < sizes.size(); i++) {
        if (pt == MESC && i == primary_part) continue;
        if (sizes[i] > limit - sum) {
            throw overflow_error("Star-clique objective sum exceeds uint64");
        }
        sum += sizes[i];
    }
    if (pt == MVSC) return sum;
    if (sizes[primary_part] != 0 && sum > limit / sizes[primary_part]) {
        throw overflow_error("Star-clique MESC objective exceeds uint64");
    }
    return static_cast<ull>(sizes[primary_part]) * sum;
}

#endif
