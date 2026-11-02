#include "pattern_core.h"
#include "common/text.h"
#include <algorithm>

namespace mcviz {
void CoreWorkspace::bind(const Graph& g, const Sets& domains, CoreMetrics* metrics) {
    offsets_.resize(g.K + 1);
    offsets_[0] = 0;
    for (Layer i = 0; i < g.K; ++i) offsets_[i + 1] = offsets_[i] + g.vertex_ids[i].size();
    if (stamps_.size() != offsets_.back()) {
        stamps_.assign(offsets_.back(), 0);
        ids_.resize(offsets_.back());

        epoch_ = 0;
        if (metrics) ++metrics->map_initializations;
    }
    if (++epoch_ == 0) {
        std::fill(stamps_.begin(), stamps_.end(), 0);
        epoch_ = 1;
        if (metrics) ++metrics->map_initializations;
    }
    std::size_t id = 0;
    for (Layer i = 0; i < g.K; ++i) for (auto u : domains[i]) {
        auto pos = offsets_[i] + u;
        stamps_[pos] = epoch_;
        ids_[pos] = id++;
    }
}

std::optional<std::size_t> CoreWorkspace::find(Layer i, VertexPos u) const {
    auto pos = offsets_[i] + u;

    if (stamps_[pos] != epoch_) return std::nullopt;
    return ids_[pos];
}

// Peel vertices missing a required layer degree; every feasible tuple survives.
// Peeling scans each directed working edge at most once. Counters use
// O(sum_i |C_i| deg_P(i)) space, in addition to O(n+K+|P|+m_P) mappings and arcs.
CoreResult reduce_queue(const Graph& g, const LayerPattern& p, const Thresholds& gamma,
                        Sets domains, const Required& required, CoreWorkspace* workspace,
                        CoreMetrics* metrics) {
    if (metrics) *metrics = {};

    auto conditioned = condition_domains(g, p, std::move(domains), required);
    CoreResult out{std::move(conditioned.domains), conditioned.possible, 0, 0};
    if (!out.possible) return out;

    Thresholds lower = gamma;

    Sets mandatory(g.K);
    for (auto q : required) mandatory[q.layer].push_back(q.pos);

    for (Layer i = 0; i < g.K; ++i) {
        auto& a = mandatory[i];

        std::sort(a.begin(), a.end());
        a.erase(std::unique(a.begin(), a.end()), a.end());
        lower[i] = std::max<Score>(lower[i], a.size());
    }

    CoreWorkspace local;
    auto& map = workspace ? *workspace : local;
    map.bind(g, out.domains, metrics);
    std::vector<LocalVertex> vertices;

    std::vector<std::size_t> base{0};
    for (Layer i = 0; i < g.K; ++i) for (auto u : out.domains[i]) {
        vertices.push_back({i, u});
        base.push_back(base.back() + p.neighbors[i].size());
    }
    const auto n = vertices.size();
    std::vector<Score> deg(base.back(), 0);
    struct Arc { std::size_t v, slot; };
    std::vector<Arc> arcs;

    std::vector<std::size_t> start(n + 1), bad(n, 0);

    std::vector<std::vector<std::size_t>> reverse(g.K);
    for (Layer i = 0; i < g.K; ++i) for (auto j : p.neighbors[i]) {
        auto it = std::lower_bound(p.neighbors[j].begin(), p.neighbors[j].end(), i);
        reverse[i].push_back(it - p.neighbors[j].begin());
    }

    for (std::size_t u = 0; u < n; ++u) {
        auto x = vertices[u];
        start[u] = arcs.size();
        for (std::size_t s = 0; s < p.neighbors[x.layer].size(); ++s) {
            auto j = p.neighbors[x.layer][s];
            for (auto v : g.neighbors(x.layer, x.pos, j)) {
                if (metrics) ++metrics->adjacency_scans;
                auto id = map.find(j, v);
                if (!id) continue;

                ++deg[base[u] + s];

                arcs.push_back({*id, base[*id] + reverse[x.layer][s]});
            }

            if (deg[base[u] + s] < lower[j]) ++bad[u];
        }
    }
    start[n] = arcs.size();
    if (metrics) { metrics->slots = deg.size(); metrics->arcs = arcs.size(); }

    std::vector<unsigned char> alive(n, 1), queued(n, 0);
    std::vector<std::size_t> queue;
    for (std::size_t u = 0; u < n; ++u) if (bad[u]) {
        queue.push_back(u);
        queued[u] = 1;
    }

    for (std::size_t head = 0; head < queue.size(); ++head) {
        auto u = queue[head];

        alive[u] = 0;
        ++out.removed;

        if (metrics) ++metrics->queue_pops;
        auto threshold = lower[vertices[u].layer];
        for (auto e = start[u]; e < start[u + 1]; ++e) {
            auto a = arcs[e];

            if (!alive[a.v]) continue;
            auto old = deg[a.slot]--;
            if (metrics) ++metrics->degree_updates;
            if (old == threshold) {
                ++bad[a.v];

                if (!queued[a.v]) {
                    queued[a.v] = 1;
                    queue.push_back(a.v);
                }
            }
        }
    }
    Sets keep(g.K);
    for (std::size_t u = 0; u < n; ++u) if (alive[u]) {
        auto x = vertices[u];
        keep[x.layer].push_back(x.pos);
    }
    out.domains = std::move(keep);
    out.possible = enough_capacity(out.domains, lower) && contains_required(out.domains, required);
    return out;
}
}
