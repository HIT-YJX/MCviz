#include "api.h"
#include "bcvizm_lp.h"
#include "starclique.h"
#include <algorithm>

namespace mcviz::baseline {
namespace {
ProblemType legacy_objective(Objective obj) {
    if (obj == Objective::edge) return MESC;
    if (obj == Objective::vertex) return MVSC;
    return MBSC;
}

std::vector<std::vector<bool>> matrix(const LayerPattern& pattern) {
    std::vector<std::vector<bool>> P(pattern.K, std::vector<bool>(pattern.K, false));
    for (auto [i, j] : pattern.edges) P[i][j] = P[j][i] = true;
    return P;
}

}

SearchResult run(const Graph& graph, const LayerPattern& pattern, Objective objective,
                 const QueryRequest& request, const std::string& method,
                 unsigned) {
    if (method == "bcvizm") return run_bcvizm(graph, pattern, objective, request);

    SearchResult out;
    out.sets.resize(graph.K);
    std::optional<LocalVertex> q;
    if (request.required) {
        q = LocalVertex{request.required->layer,
                        graph.get_pos(request.required->layer, request.required->id)};
    }
    std::vector<uint> tau(graph.K);
    for (Layer p = 0; p < graph.K; p++) {
        if (request.gamma[p] > graph.vertex_ids[p].size()) {
            out.status = Status::infeasible;
            return out;
        }
        tau[p] = static_cast<uint>(request.gamma[p]);
    }

    KPartiteGraph local;
    local.load(graph);
    unsigned anchor = pattern.anchor_layer;
    if (anchor >= graph.K) anchor = 0;
    StarCliqueFinder finder(local, legacy_objective(objective), tau, anchor);
    finder.required = q;
    const auto P = matrix(pattern);

    StarClique result;
    if (method == "motif") result = finder.find_max_pattern_clique_motif_online(P);
    else if (method == "kpc") result = finder.find_max_pattern_clique_kpc_online(P);
    else {
        auto centers = pattern.star_centers();
        if (!q && !centers.empty()) {
            if (std::find(centers.begin(), centers.end(), anchor) == centers.end()) {
                anchor = centers.front();
                StarCliqueFinder star_finder(local, legacy_objective(objective), tau, anchor);
                result = star_finder.find_max_star_clique_kspc_online(graph.K);
            } else {
                result = finder.find_max_star_clique_kspc_online(graph.K);
            }
        } else {
            result = finder.find_max_pattern_clique_kspc_online(P, graph.K);
        }
    }

    if (!result.found) {
        out.status = Status::infeasible;
        return out;
    }
    out.found = true;
    out.score = result.size;
    out.sets = std::move(result.sets);
    out.status = Status::optimal;
    return out;
}
}
