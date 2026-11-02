#include "build_filter.h"
#include "search/pattern_filter.h"
#include "common/text.h"
#include <chrono>
#include <new>
#include <limits>
#include <stdexcept>

namespace mcviz {
namespace {
using Clock = std::chrono::steady_clock;
double ms(Clock::time_point start) {
    return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}


Sets peel_filtered(const Graph& g, const BuildSpec& spec, const Sets& core,
                   const std::vector<unsigned char>& keep) {
    std::vector<std::vector<unsigned char>> member(g.K);
    std::vector<VertexRef> vertices;
    for (Layer i = 0; i < g.K; ++i) {
        member[i].resize(g.vertex_ids[i].size());
        for (auto u : core[i]) { member[i][u] = 1; vertices.push_back({i, u}); }
    }
    std::vector<std::pair<VertexRef, VertexRef>> edges;
    for (std::size_t id = 0; id < g.edges.size(); ++id) {
        const auto& e = g.edges[id];
        if (keep[id] && spec.pattern.adjacent(e.u.layer, e.v.layer) &&
            member[e.u.layer][e.u.pos] && member[e.v.layer][e.v.pos])
            edges.push_back({{e.u.layer, e.u.pos}, {e.v.layer, e.v.pos}});
    }
    auto active = Graph::from_records(g.K, vertices, edges);
    auto reduced = reduce(active, spec.pattern, spec.beta, active.all_vertices(), {}, CoreBackend::queue);
    Sets result(g.K);
    for (Layer i = 0; i < g.K; ++i) for (auto u : reduced.domains[i])
        result[i].push_back(static_cast<VertexPos>(active.get_id(i, u)));
    return result;
}
}


BuildFilter prepare_build_filter(PatternIndex& idx, const Graph& g, const BuildSpec& spec, bool possible) {
    BuildFilter out{idx.core_domains, {}, possible};
    auto& c = idx.build_counters;
    for (auto key : {"build_filter_attempted", "build_filter_removed_edges", "build_filter_removed_vertices",
                     "build_filter_zeroed_weights", "build_filter_recore_estimated_bytes", "build_filter_recore_skipped_memory", "build_filter_recore_complete",
                     "build_filter_eligible_pairs", "build_filter_completed_pairs", "build_filter_skipped_work",
                     "build_filter_skipped_memory", "build_filter_skipped_threshold", "build_filter_work",
                     "build_filter_wedges", "build_filter_peak_bytes", "build_filter_mask_available"}) c[key] = 0;
    for (auto key : {"build_filter_ms", "build_filter_prepare_ms", "build_filter_pair_ms", "build_filter_recore_ms"})
        idx.build_timing[key] = 0;
    for (Layer i = 0; i < g.K; ++i) c["build_filter_vertices_layer_" + std::to_string(i)] = out.domains[i].size();
    if (!possible || idx.options.build_filter == "off") return out;
    const auto start = Clock::now();
    c["build_filter_attempted"] = 1;
    auto result = filter_butterflies(g, spec.pattern, spec.beta, out.domains,
        idx.options.build_filter_work_budget, idx.options.build_filter_cache_bytes);
    c["build_filter_removed_edges"] = result.removed_edges;
    c["build_filter_eligible_pairs"] = result.eligible_pairs;
    c["build_filter_completed_pairs"] = result.completed_pairs;
    c["build_filter_skipped_work"] = result.skipped_work;
    c["build_filter_skipped_memory"] = result.skipped_memory;
    c["build_filter_skipped_threshold"] = result.skipped_threshold;
    c["build_filter_work"] = result.work; c["build_filter_wedges"] = result.wedges;
    c["build_filter_peak_bytes"] = result.peak_bytes;
    c["build_filter_mask_available"] = result.mask_available;
    idx.build_timing["build_filter_prepare_ms"] = result.prepare_ms;
    idx.build_timing["build_filter_pair_ms"] = result.filter_ms;
    out.keep_edges = std::move(result.keep_edges);
    if (result.removed_edges) {
        std::uint64_t n = 0, original = 0;
        for (Layer i = 0; i < g.K; ++i) { n += out.domains[i].size(); original += g.vertex_ids[i].size(); }

        long double estimate = out.keep_edges.size() + original + 512.0L * (n + g.edges.size() + spec.pattern.edges.size() + 1);
        c["build_filter_recore_estimated_bytes"] = estimate >= std::numeric_limits<std::uint64_t>::max()
            ? std::numeric_limits<std::uint64_t>::max() : static_cast<std::uint64_t>(estimate);
        if (estimate > idx.options.build_filter_cache_bytes) c["build_filter_recore_skipped_memory"] = 1;
        else {
            const auto recore = Clock::now();
            {
                auto domains = peel_filtered(g, spec, out.domains, out.keep_edges);
                std::uint64_t left = 0;
                for (const auto& d : domains) left += d.size();
                c["build_filter_removed_vertices"] = n - left;
                out.domains = std::move(domains);
                out.possible = enough_capacity(out.domains, spec.beta);
                c["build_filter_recore_complete"] = 1;
            }
            idx.build_timing["build_filter_recore_ms"] = ms(recore);
        }
    }
    for (Layer i = 0; i < g.K; ++i) c["build_filter_vertices_layer_" + std::to_string(i)] = out.domains[i].size();
    idx.build_timing["build_filter_ms"] = ms(start);
    return out;
}


void apply_build_filter(PatternIndex& idx, const Graph& g, const BuildSpec& spec, const BuildFilter& filter) {
    if (filter.keep_edges.empty()) return;
    const auto started = Clock::now();
    std::size_t at = 0;
    for (std::size_t id = 0; id < g.edges.size(); ++id) {
        const auto& e = g.edges[id];
        if (!spec.pattern.adjacent(e.u.layer, e.v.layer)) continue;
        if (!filter.keep_edges[id]) {
            if (idx.weights[at].weight) ++idx.build_counters["build_filter_zeroed_weights"];
            idx.weights[at].weight = 0;
        }
        ++at;
    }
    idx.build_timing["build_filter_ms"] += ms(started);
}
}
