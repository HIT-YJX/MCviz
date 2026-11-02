#include "pattern_core.h"
#include "common/text.h"
#include <algorithm>

namespace mcviz {
CoreResult reduce(const Graph& g, const LayerPattern& p, const Thresholds& gamma,
                  Sets domains, const Required& required, CoreBackend backend,
                  CoreWorkspace* workspace, CoreMetrics* metrics) {
    if (metrics) *metrics = {};
    if (backend == CoreBackend::queue)
        return reduce_queue(g, p, gamma, std::move(domains), required, workspace, metrics);

    auto conditioned = condition_domains(g, p, std::move(domains), required);

    Thresholds lower = gamma;
    Sets mandatory(g.K);

    for (auto q : required) mandatory[q.layer].push_back(q.pos);
    for (Layer i = 0; i < g.K; ++i) {
        auto& a = mandatory[i];
        std::sort(a.begin(), a.end());
        a.erase(std::unique(a.begin(), a.end()), a.end());
        lower[i] = std::max<Score>(lower[i], a.size());
    }
    CoreResult out{std::move(conditioned.domains), conditioned.possible, 0, 0};
    if (!out.possible) return out;
    while (true) {
        Sets keep(g.K);
        std::uint64_t removed = 0;

        for (Layer i = 0; i < g.K; ++i) for (auto u : out.domains[i]) {
            bool valid = true;
            for (auto j : p.neighbors[i]) {
                Score deg = 0;
                for (auto v : out.domains[j]) {
                    if (metrics) ++metrics->pair_probes;
                    if (g.has_edge(i, u, j, v)) ++deg;
                }
                if (deg < lower[j]) { valid = false; break; }
            }

            if (valid) keep[i].push_back(u);
            else ++removed;
        }
        ++out.rounds;
        out.removed += removed;

        out.domains = std::move(keep);
        out.possible = enough_capacity(out.domains, lower) && contains_required(out.domains, required);
        if (!out.possible || removed == 0) return out;
    }
}
}
