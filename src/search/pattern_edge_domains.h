#pragma once
#include "pattern_core.h"
#include <memory>

namespace mcviz {
class PatternAdjacency;

struct PatternEdgeDomainCounters {
    std::uint64_t requests = 0, cache_requests = 0, core_calls = 0, impossible = 0;
    std::uint64_t sparse_requests = 0;
    std::uint64_t base_vertices = 0, condition_vertices = 0, output_vertices = 0;
    std::uint64_t condition_removed = 0, core_removed = 0;
    std::uint64_t adjacency_scans = 0, degree_updates = 0, queue_pops = 0;
    std::uint64_t map_initializations = 0, peak_domain_vertices = 0;
    std::uint64_t peak_core_slots = 0, peak_core_arcs = 0;
    // These counts include only storage owned by this workspace.
    std::uint64_t adjacency_logical_bytes = 0, persistent_logical_bytes = 0;
    // Informational only: externally shared cache payload is charged once by
    // its owner, and excluded from this workspace's persistent/peak estimates.
    std::uint64_t shared_adjacency_logical_bytes = 0;
    // Queue internals expose element counts, not allocated capacities. These are
    // conservative element-storage estimates, excluding allocator slack and Graph.
    std::uint64_t workspace_bytes_estimate = 0, core_peak_bytes_estimate = 0;
    std::uint64_t peak_bytes_estimate = 0;
};

struct PatternEdgeDomainTiming {
    double prepare_ms = 0, cache_prepare_ms = 0;
    double condition_ms = 0, core_ms = 0, get_ms = 0;
};

// A serial batch workspace for arbitrary normalized P and positive beta. It owns
// beta/base copies; Graph and P must outlive it and remain unchanged. No star,
// connectedness or bipartiteness assumption is made. A shared adjacency also
// borrows its original base; its owner must keep that input alive until all
// workspaces finish. Each thread needs a separate PatternEdgeDomains instance.
class PatternEdgeDomains {
public:
    // With shared adjacency, cache/cache_bytes are unused and cache_prepare_ms
    // is zero; that object's owner applies its preparation budget/accounting.
    PatternEdgeDomains(const Graph& graph, const LayerPattern& pattern,
                       const Thresholds& beta, const Sets& base,
                       bool cache = true, bool local_core = true,
                       std::uint64_t cache_bytes = 64ULL * 1024 * 1024,
                       std::shared_ptr<const PatternAdjacency> shared = {});
    ~PatternEdgeDomains();
    PatternEdgeDomains(const PatternEdgeDomains&) = delete;
    PatternEdgeDomains& operator=(const PatternEdgeDomains&) = delete;

    // Condition all K domains on the two endpoints; optionally peel their queue
    // beta-core with both endpoints mandatory. possible=false proves kappa(e)=0;
    // possible=true is a necessary condition, never a positive-cohesion certificate.
    // The complete returned domains preserve every beta-witness containing e.
    CoreResult get(const GraphEdge& edge);
    bool cache_ready() const;
    const PatternEdgeDomainCounters& counters() const { return counters_; }
    const PatternEdgeDomainTiming& timing() const { return timing_; }

private:
    const Graph& graph_;
    const LayerPattern& pattern_;
    Thresholds beta_;
    Sets base_;
    bool local_core_;
    // When supplied externally, cache/cache_bytes do not rebuild or change it.
    // Its owner accounts preparation time and storage once across all workers.
    std::shared_ptr<const PatternAdjacency> adjacency_;
    CoreWorkspace workspace_;
    PatternEdgeDomainCounters counters_;
    PatternEdgeDomainTiming timing_;
    std::uint64_t original_vertices_ = 0;
};
} // namespace mcviz
