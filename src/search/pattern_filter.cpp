#include "pattern_filter.h"
#include "common/text.h"
#include <algorithm>
#include <chrono>
#include <limits>
#include <new>
#include <stdexcept>

namespace mcviz {
namespace {
using Clock = std::chrono::steady_clock;
enum class Refusal { none, work, memory };

std::uint64_t sum(std::uint64_t a, std::uint64_t b) {
    return b > std::numeric_limits<std::uint64_t>::max() - a
        ? std::numeric_limits<std::uint64_t>::max() : a + b;
}

class PairCounts {
    std::vector<std::uint64_t> keys;
    std::vector<unsigned char> counts;
    std::size_t used = 0;
public:
    void allocate(std::size_t n) { keys.assign(n, 0); counts.assign(n, 0); }

    template<class Charge>
    unsigned count(std::uint64_t key, bool insert, Charge&& charge) {
        auto h = key + 0x9e3779b97f4a7c15ULL;
        h = (h ^ (h >> 30)) * 0xbf58476d1ce4e5b9ULL;
        h = (h ^ (h >> 27)) * 0x94d049bb133111ebULL;

        auto at = std::size_t(h ^ (h >> 31)) & (keys.size() - 1);
        for (;;) {
            if (!charge(1)) return 3;

            if (!keys[at]) {
                if (!insert) return 0;
                if (used >= keys.size() / 2) return 4;
                keys[at] = key; counts[at] = 1; ++used;
                return 1;
            }
            if (keys[at] == key) {
                if (insert && counts[at] < 2) ++counts[at];
                return counts[at];
            }
            at = (at + 1) & (keys.size() - 1);
        }
    }
};

class Filter {
    const Graph& g;
    const Sets& core;
    PatternFilterResult& out;
    std::uint64_t budget, memory, base_bytes = 0;
    std::vector<std::vector<unsigned char>> member;
    std::vector<std::size_t> edge_start;

    bool charge(std::uint64_t amount) {
        if (amount > budget - out.work) { refusal = Refusal::work; return false; }
        out.work += amount; return true;
    }



    // Track estimated memory and refuse work exceeding the configured budget.
    bool check_memory(long double bytes) {
        auto value = bytes >= std::numeric_limits<std::uint64_t>::max()
            ? std::numeric_limits<std::uint64_t>::max() : static_cast<std::uint64_t>(bytes);

        out.peak_bytes = std::max(out.peak_bytes, value);

        if (bytes > memory) { refusal = Refusal::memory; return false; }

        return true;
    }

    bool live(const GraphEdge& e, Layer j) const {
        return e.v.layer == j && member[e.u.layer][e.u.pos] && member[j][e.v.pos];
    }


    VertexPos middle(std::size_t id, bool reverse) const {
        const auto& e = g.edges[id];
        return reverse ? e.v.pos : e.u.pos;
    }



    VertexPos other(std::size_t id, bool reverse) const {
        const auto& e = g.edges[id];
        return reverse ? e.u.pos : e.v.pos;
    }
public:
    Refusal refusal = Refusal::none;
    Filter(const Graph& graph, const Sets& domains, PatternFilterResult& result,
           std::uint64_t work, std::uint64_t bytes)
        : g(graph), core(domains), out(result), budget(work), memory(bytes) {}

    void prepare() {
        if (!charge(1)) return;
        std::uint64_t n = 0;
        for (const auto& v : g.vertex_ids) n = sum(n, v.size());

        base_bytes = sum(sum(n, g.edges.size()), 64ULL * (g.K + 1));

        if (!check_memory(base_bytes)) return;
        if (!charge(sum(n, g.edges.size()))) return;
        member.resize(g.K);

        for (Layer i = 0; i < g.K; ++i) {
            member[i].assign(g.vertex_ids[i].size(), 0);
            for (auto u : core[i]) {
                member[i][u] = 1;
            }
        }

        edge_start.assign(g.K + 1, 0);
        for (const auto& e : g.edges) ++edge_start[e.u.layer + 1];

        for (Layer i = 0; i < g.K; ++i) edge_start[i + 1] += edge_start[i];
        out.keep_edges.assign(g.edges.size(), 1);
        out.mask_available = true;
    }



    std::pair<std::uint64_t, std::uint64_t> volume(Layer i, Layer j) {
        std::uint64_t m = 0, w = 0;
        if (!charge(core[i].size())) return {};
        for (auto u : core[i]) {
            const auto& adj = g.neighbors(i, u, j);
            if (!charge(adj.size())) return {};
            std::uint64_t d = 0;
            for (auto v : adj) d += member[j][v];

            m = sum(m, d); w = sum(w, d * (d ? d - 1 : 0) / 2);
        }
        return {m, w};
    }


    void pair(Layer i, Layer j) {
        refusal = Refusal::none;
        const auto left = volume(i, j); if (refusal != Refusal::none) return;

        const auto right = volume(j, i); if (refusal != Refusal::none) return;

        const auto m = left.first, w = std::min(left.second, right.second);
        out.wedges = sum(out.wedges, w);
        if (!m) { ++out.completed_pairs; return; }

        if (w > (budget - out.work) / 2) { refusal = Refusal::work; return; }

        const bool reverse = right.second < left.second;
        const auto rows = g.vertex_ids[reverse ? j : i].size();
        if (!check_memory(static_cast<long double>(base_bytes) + static_cast<long double>(sizeof(std::size_t)) * (m + rows + 1))) return;
        if (!charge(sum(rows + 1, 2 * (edge_start[i + 1] - edge_start[i])))) return;
        std::vector<std::size_t> offsets(rows + 1, 0), ids(m);
        for (auto e = edge_start[i]; e < edge_start[i + 1]; ++e) if (live(g.edges[e], j))
            ++offsets[middle(e, reverse) + 1];

        for (std::size_t u = 0; u < rows; ++u) offsets[u + 1] += offsets[u];
        for (auto e = edge_start[i]; e < edge_start[i + 1]; ++e) if (live(g.edges[e], j))
            ids[--offsets[middle(e, reverse) + 1]] = e;
        std::vector<std::size_t>().swap(offsets);

        const long double fixed = static_cast<long double>(base_bytes) + (sizeof(std::size_t) + 1.0L) * m + 256;
        if (!check_memory(fixed)) return;

        std::vector<unsigned char> supported(m, 0);
        if (w) {
            const auto available = static_cast<std::uint64_t>((memory - fixed) / 9);
            std::size_t slots = 1;
            const auto desired = sum(w, w);
            while (slots <= available / 2 && slots < desired) slots *= 2;
            if (slots < 2) { refusal = Refusal::memory; return; }

            if (!check_memory(fixed + 9.0L * slots)) return;
            PairCounts counts; counts.allocate(slots);
            auto debit = [&](std::uint64_t n) { return charge(n); };

            for (unsigned pass = 0; pass < 2; ++pass) {
                if (!charge(m)) return;
                for (std::size_t first = 0; first < ids.size();) {
                    auto last = first + 1;
                    while (last < ids.size() && middle(ids[last], reverse) == middle(ids[first], reverse)) ++last;
                    for (auto a = first; a < last; ++a) for (auto b = a + 1; b < last; ++b) {
                        auto x = other(ids[a], reverse), y = other(ids[b], reverse);
                        if (x > y) std::swap(x, y);
                        const auto key = (std::uint64_t(x) << 32) | y;
                        auto count = counts.count(key, pass == 0, debit);
                        if (count == 3) return;

                        if (count == 4) { refusal = Refusal::memory; return; }
                        if (pass && count == 2) supported[a] = supported[b] = 1;
                    }
                    first = last;
                }
            }
        }

        if (!charge(m)) return;

        for (std::size_t e = 0; e < ids.size(); ++e) if (!supported[e]) {
            out.keep_edges[ids[e]] = 0;
            ++out.removed_edges;
        }
        ++out.completed_pairs;
    }
};
}

// Butterfly deficits safely reject pattern edges when both layer thresholds are at least two.
// Work accounting and cache admission follow their budgets; refusal preserves unprocessed edges.
PatternFilterResult filter_butterflies(const Graph& g, const LayerPattern& p,
                                      const Thresholds& beta, const Sets& core,
                                      std::uint64_t work_budget, std::uint64_t memory_budget) {
    PatternFilterResult out;

    for (auto [i, j] : p.edges) {
        if (beta[i] >= 2 && beta[j] >= 2) ++out.eligible_pairs;
        else ++out.skipped_threshold;
    }
    if (!out.eligible_pairs) return out;
    auto started = Clock::now();
    Filter trial(g, core, out, work_budget, memory_budget);
    trial.prepare();
    bool prepared = trial.refusal == Refusal::none;
    if (trial.refusal == Refusal::work) out.skipped_work = out.eligible_pairs;
    if (trial.refusal == Refusal::memory) out.skipped_memory = out.eligible_pairs;

    out.prepare_ms = std::chrono::duration<double, std::milli>(Clock::now() - started).count();
    if (!prepared) return out;
    started = Clock::now();
    for (auto [i, j] : p.edges) if (beta[i] >= 2 && beta[j] >= 2) {
        trial.pair(i, j);

        if (trial.refusal == Refusal::work) ++out.skipped_work;
        if (trial.refusal == Refusal::memory) ++out.skipped_memory;
    }
    out.filter_ms = std::chrono::duration<double, std::milli>(Clock::now() - started).count();
    return out;
}
}
