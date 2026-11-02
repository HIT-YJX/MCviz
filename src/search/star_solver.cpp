#include "star_solver.h"
#include "common/text.h"
#include <algorithm>
#include <iterator>

namespace mcviz {
// Choose the requested anchor when it is a valid star center.
static Layer center_of(const LayerPattern& p) {
    auto centers = p.star_centers();

    return std::find(centers.begin(), centers.end(), p.anchor_layer) == centers.end() ? centers[0] : p.anchor_layer;
}

namespace {
struct ReferenceArithmetic {
    // Add two scores through the shared overflow-aware arithmetic helper.
    Score add(Score a, Score b) const { return checked_add(a, b); }


    Score multiply(Score a, Score b) const { return checked_mul(a, b); }


    Score capacity(const LayerPattern& p, Objective obj, const std::vector<Score>& cap) const {
        return score(p, obj, cap);
    }
};

struct PreparedArithmetic {
    Score add(Score a, Score b) const {
        if (b > std::numeric_limits<Score>::max() - a) throw std::runtime_error("Score addition overflow");
        return a + b;
    }



    Score multiply(Score a, Score b) const {
        if (a != 0 && b > std::numeric_limits<Score>::max() / a) throw std::runtime_error("Score multiplication overflow");
        return a * b;
    }


    // Evaluate the objective bound induced by the supplied layer capacities.
    Score capacity(const LayerPattern& p, Objective obj, const std::vector<Score>& cap) const {
        Score value = 0;
        if (obj == Objective::edge) {
            for (auto edge : p.edges) value = add(value, multiply(cap[edge.first], cap[edge.second]));
        } else if (obj == Objective::vertex) {
            for (auto size : cap) value = add(value, size);
        } else value = *std::min_element(cap.begin(), cap.end());
        return value;
    }
};

template<class Adjacent, class Arithmetic>
Score evaluate_star_bound(unsigned layers, const LayerPattern& p, Objective obj,
                          const Thresholds& beta, GraphEdge e, bool support, Layer center,
                          Adjacent adjacency, Arithmetic arithmetic) {
    const auto& A = adjacency(e.v.layer, e.v.pos, center);
    std::vector<Score> cap(layers);
    cap[center] = A.size();
    for (Layer j = 0; j < layers; ++j) if (j != center) cap[j] = adjacency(center, e.u.pos, j).size();
    for (Layer j = 0; j < layers; ++j) if (cap[j] < beta[j]) return 0;
    Score coarse = arithmetic.capacity(p, obj, cap);
    if (!support) return coarse;

    std::vector<Score> sums(A.size() + 1, 0), minima(A.size() + 1, std::numeric_limits<Score>::max());
    std::vector<bool> feasible(A.size() + 1, true);
    for (Layer j = 0; j < layers; ++j) if (j != center) {
        std::vector<Score> hist(A.size() + 1, 0);
        for (auto v : adjacency(center, e.u.pos, j)) {
            const auto& B = adjacency(j, v, center);
            std::size_t ai = 0, bi = 0, count = 0;
            while (ai < A.size() && bi < B.size()) {
                if (A[ai] < B[bi]) ++ai;
                else if (B[bi] < A[ai]) ++bi;
                else { if (A[ai] != e.u.pos) ++count; ++ai; ++bi; }
            }
            ++hist[count];
        }
        Score count = 0;
        for (std::size_t a = A.size(); a > 0; --a) {
            count += hist[a - 1]; sums[a] = arithmetic.add(sums[a], count);
            minima[a] = std::min(minima[a], count);

            if (count < beta[j]) feasible[a] = false;
        }
    }
    Score bound = 0;
    for (std::size_t a = static_cast<std::size_t>(beta[center]); a <= A.size(); ++a) if (feasible[a]) {
        Score value = obj == Objective::edge ? arithmetic.multiply(a, sums[a]) :
            (obj == Objective::vertex ? arithmetic.add(a, sums[a]) : std::min<Score>(a, minima[a]));
        bound = std::max(bound, value);
    }
    return std::min(coarse, bound);
}
}

Score star_edge_bound(const Graph& g, const LayerPattern& p, Objective obj,
                      const Thresholds& beta, const Sets& core, GraphEdge e, bool support) {
    Layer center = center_of(p);
    if (e.u.layer != center) std::swap(e.u, e.v);

    auto adjacency = [&](Layer i, VertexPos u, Layer j) {
        std::vector<VertexPos> a;
        const auto& ns = g.neighbors(i, u, j);

        std::set_intersection(ns.begin(), ns.end(), core[j].begin(), core[j].end(), std::back_inserter(a));
        return a;
    };
    return evaluate_star_bound(g.K, p, obj, beta, e, support, center, adjacency, ReferenceArithmetic{});
}



PreparedStarBounds::PreparedStarBounds(const Graph& active, const LayerPattern& pattern,
                                     Objective objective, const Thresholds& beta)
    : pattern_(pattern), beta_(beta), objective_(objective), center_(center_of(pattern)),
      layers_(active.K), primary_count_(active.vertex_ids.at(center_).size()) {
    require(primary_count_ <= std::numeric_limits<std::size_t>::max() / layers_, "Prepared star adjacency size overflow");
    primary_adjacency_.resize(primary_count_ * layers_, nullptr);
    secondary_adjacency_.resize(layers_);
    for (VertexPos u = 0; u < primary_count_; ++u) for (Layer j = 0; j < layers_; ++j) if (j != center_) {
        const auto& neighbors = active.neighbors(center_, u, j);

        primary_adjacency_[static_cast<std::size_t>(u) * layers_ + j] = &neighbors;
    }
    for (Layer j = 0; j < layers_; ++j) if (j != center_) {
        auto& slots = secondary_adjacency_[j];
        slots.resize(active.vertex_ids[j].size());
        for (VertexPos v = 0; v < slots.size(); ++v) {
            const auto& neighbors = active.neighbors(j, v, center_);

            slots[v] = &neighbors;
        }
    }
}


const PreparedStarBounds::Neighborhood& PreparedStarBounds::adjacency(Layer layer, VertexPos vertex,
                                                                     Layer target) const {
    return layer == center_ ? *primary_adjacency_[static_cast<std::size_t>(vertex) * layers_ + target]
                           : *secondary_adjacency_[layer][vertex];
}

Score PreparedStarBounds::bound(GraphEdge edge, bool support) const {
    if (edge.u.layer != center_) std::swap(edge.u, edge.v);
    auto neighbors = [&](Layer i, VertexPos u, Layer j) -> const Neighborhood& { return adjacency(i, u, j); };
    return evaluate_star_bound(layers_, pattern_, objective_, beta_, edge, support, center_, neighbors, PreparedArithmetic{});
}
}
