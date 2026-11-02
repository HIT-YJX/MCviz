#include "pattern_edge_domains.h"
#include "pattern_adjacency.h"
#include "common/text.h"
#include <algorithm>
#include <chrono>
#include <limits>

namespace mcviz {
namespace {
using Clock = std::chrono::steady_clock;

// Account successful and exceptional calls without changing their result; O(1).
struct Timer {
    double& value;
    Clock::time_point start = Clock::now();
    explicit Timer(double& ms) : value(ms) {}


    ~Timer() { value += std::chrono::duration<double, std::milli>(Clock::now() - start).count(); }
};

// Saturation applies only to statistics/space estimates, never objective values.
std::uint64_t add(std::uint64_t a, std::uint64_t b) {
    const auto limit = std::numeric_limits<std::uint64_t>::max();
    return b > limit - a ? limit : a + b;
}

// Counted arrays may exceed the reporting range; O(1), no allocation.
std::uint64_t mul(std::uint64_t n, std::uint64_t bytes) {
    const auto limit = std::numeric_limits<std::uint64_t>::max();

    return bytes && n > limit / bytes ? limit : n * bytes;
}

// Domain cardinality includes every layer, including isolated pattern layers.
// O(K) time, O(1) scratch.
std::uint64_t vertices(const Sets& d) {
    std::uint64_t n = 0;
    for (const auto& a : d) n = add(n, a.size());
    return n;
}

// Conservative queue element storage: two domain copies, typed vertices, CSR,
// degrees, queue/flags, mandatory sets and reverse pattern slots. Vector capacity
// slack and allocator metadata are unavailable through CoreMetrics and excluded.
// O(1) time/space; n is the conditioned size, not the original Graph size.
std::uint64_t core_bytes(unsigned k, std::size_t edges, std::uint64_t n,
                         const CoreMetrics& metrics) {
    auto bytes = mul(4ULL * k, sizeof(std::vector<VertexPos>));
    bytes = add(bytes, mul(k, sizeof(Score)));
    bytes = add(bytes, mul(n, sizeof(LocalVertex) + 2 * sizeof(VertexPos) + 2));
    bytes = add(bytes, mul(add(mul(4, n), 2), sizeof(std::size_t)));
    bytes = add(bytes, mul(metrics.slots, sizeof(Score)));

    bytes = add(bytes, mul(metrics.arcs, 2 * sizeof(std::size_t)));
    bytes = add(bytes, mul(mul(2, edges), sizeof(std::size_t)));

    return add(bytes, 2 * (sizeof(VertexPos) + sizeof(LocalVertex)));
}
} // namespace

// Normalize/own the base once, then prepare or borrow typed singleton adjacency.
// Beyond cache preparation: O(sum_i |C_i|log(|C_i|+1)+K) time, O(n_C+K) storage.
// External storage/time is charged by its owner, while queue state stays private.
PatternEdgeDomains::PatternEdgeDomains(const Graph& g, const LayerPattern& p,
                                     const Thresholds& beta, const Sets& base,
                                     bool cache, bool local_core, std::uint64_t cache_bytes,
                                     std::shared_ptr<const PatternAdjacency> shared)
    : graph_(g), pattern_(p), local_core_(local_core) {
    Timer timer(timing_.prepare_ms);
    beta_ = beta;
    base_ = condition_domains(g, p, base, {}).domains;
    counters_.base_vertices = vertices(base_);
    for (const auto& a : g.vertex_ids) original_vertices_ = add(original_vertices_, a.size());
    auto persistent = add(sizeof(*this), mul(beta_.capacity(), sizeof(Score)));
    persistent = add(persistent, mul(base_.capacity(), sizeof(std::vector<VertexPos>)));

    for (const auto& a : base_) persistent = add(persistent, mul(a.capacity(), sizeof(VertexPos)));

    if (shared) {
        adjacency_ = std::move(shared);
        counters_.shared_adjacency_logical_bytes = adjacency_->logical_bytes();
    } else {
        adjacency_ = std::make_shared<PatternAdjacency>(g, p, base_, cache, cache_bytes);
        timing_.cache_prepare_ms = adjacency_->prepare_ms();
        counters_.adjacency_logical_bytes = adjacency_->logical_bytes();
        persistent = add(persistent, std::max<std::uint64_t>(sizeof(PatternAdjacency), adjacency_->logical_bytes()));
    }
    counters_.persistent_logical_bytes = persistent;
    counters_.peak_bytes_estimate = persistent;
}

// Destruction releases this batch's base/workspace and its cache reference.
PatternEdgeDomains::~PatternEdgeDomains() = default;

// Read-only readiness; O(1). A missing cache never makes a domain incomplete.
bool PatternEdgeDomains::cache_ready() const { return adjacency_->ready(); }


// Conditioning and queue peeling each retain every witness containing the edge.
// The first deleted witness vertex would still have all its required neighbors,
// contradicting its deficient degree. This argument uses arbitrary P, not its shape.
// Cost is one prepared/reference conditioning plus optional reduce_queue, including
// its repeated normalization/conditioning. Workspace map initialization is amortized
// across this serial batch; every returned domain remains explicit and complete.
CoreResult PatternEdgeDomains::get(const GraphEdge& e) {
    Timer total(timing_.get_ms);

    graph_.get_id(e.u.layer, e.u.pos);
    graph_.get_id(e.v.layer, e.v.pos);

    Required required{e.u, e.v};
    counters_.requests = add(counters_.requests, 1);
    if (adjacency_->ready() && contains_required(base_, required))
        counters_.cache_requests = add(counters_.cache_requests, 1);
    else counters_.sparse_requests = add(counters_.sparse_requests, 1);

    CoreResult out;
    {
        Timer timer(timing_.condition_ms);
        auto d = adjacency_->condition(required);
        out.domains = std::move(d.domains);
        out.possible = d.possible;
    }
    const auto n = vertices(out.domains);
    counters_.condition_vertices = add(counters_.condition_vertices, n);
    counters_.condition_removed = add(counters_.condition_removed, counters_.base_vertices - n);
    counters_.peak_domain_vertices = std::max(counters_.peak_domain_vertices, n);

    // The existing queue API reconditions Q itself. Keep Q explicit so its required
    // checks and threshold contract remain exactly the production API's behavior.
    CoreMetrics metrics;
    if (local_core_) {
        counters_.core_calls = add(counters_.core_calls, 1);
        {
            Timer timer(timing_.core_ms);
            out = reduce_queue(graph_, pattern_, beta_, std::move(out.domains), required,
                               &workspace_, &metrics);
        }
        counters_.core_removed = add(counters_.core_removed, out.removed);
        counters_.adjacency_scans = add(counters_.adjacency_scans, metrics.adjacency_scans);
        counters_.degree_updates = add(counters_.degree_updates, metrics.degree_updates);

        counters_.queue_pops = add(counters_.queue_pops, metrics.queue_pops);
        counters_.map_initializations = add(counters_.map_initializations, metrics.map_initializations);
        counters_.peak_core_slots = std::max(counters_.peak_core_slots, metrics.slots);
        counters_.peak_core_arcs = std::max(counters_.peak_core_arcs, metrics.arcs);
        if (metrics.map_initializations) {
            auto bytes = mul(original_vertices_, sizeof(std::size_t) + sizeof(std::uint64_t));
            counters_.workspace_bytes_estimate = add(bytes, mul(std::uint64_t(graph_.K) + 1, sizeof(std::size_t)));
        }
        counters_.core_peak_bytes_estimate = std::max(counters_.core_peak_bytes_estimate,
            core_bytes(graph_.K, pattern_.edges.size(), n, metrics));
    }
    out.possible = out.possible && enough_capacity(out.domains, beta_)
        && contains_required(out.domains, required);
    counters_.impossible = add(counters_.impossible, !out.possible);
    counters_.output_vertices = add(counters_.output_vertices, vertices(out.domains));

    // Explicit output and a conservative condition scratch allowance are separate
    // from the adjacency cache. These estimates do not impose a process RSS limit.
    auto condition_bytes = mul(counters_.base_vertices, 2 * sizeof(VertexPos) + sizeof(std::uint64_t));

    condition_bytes = add(condition_bytes, mul(graph_.K,
        2 * sizeof(std::vector<VertexPos>) + sizeof(std::uint64_t) + sizeof(unsigned char)));
    auto temp = std::max(condition_bytes, counters_.core_peak_bytes_estimate);
    auto peak = add(counters_.persistent_logical_bytes, counters_.workspace_bytes_estimate);

    counters_.peak_bytes_estimate = std::max(counters_.peak_bytes_estimate, add(peak, temp));
    return out;
}
} // namespace mcviz
