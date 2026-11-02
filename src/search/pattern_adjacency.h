#pragma once
#include "pattern_bounds.h"

namespace mcviz {
// Reusable singleton adjacency for arbitrary P, including disconnected/empty
// patterns. The graph, normalized pattern and base must outlive this object and
// stay unchanged. Positions in all inputs/outputs are original Graph positions.
//
// Cache admission is all-or-nothing. Its logical budget includes the object,
// dynamic vector headers, sorted base copy, full original-to-base maps and bit
// arrays. Allocator metadata and per-condition output/scratch are excluded.
// Unavailable bits use sparse neighbor intersections against a base normalized
// once at construction. Unconstrained layers retain the complete base. This
// sparse representation is also valid for required points outside the base.
// After construction all methods are read-only. Concurrent condition calls own
// their output and scratch; the borrowed inputs must remain alive and unchanged.
class PatternAdjacency {
public:
    PatternAdjacency(const Graph& graph, const LayerPattern& pattern, const Sets& base,
                     bool enabled = true, std::uint64_t cache_bytes = 64ULL * 1024 * 1024);
    DomainResult condition(const Required& required) const;
    // A shared cache may only serve its original Graph/P and the same full base.
    // The supplied base is sorted by Graph position, as in condition_domains.
    bool matches(const Graph& graph, const LayerPattern& pattern, const Sets& base) const;
    // Report whether the cached adjacency representation is available for queries.
    bool ready() const { return ready_; }


    std::uint64_t logical_bytes() const { return logical_bytes_; }


    double prepare_ms() const { return prepare_ms_; }

private:
    struct Block {
        Layer from = 0, to = 0;
        std::size_t words = 0;
        std::vector<std::uint64_t> bits;
    };
    const Graph& graph_;
    const LayerPattern& pattern_;
    const Sets& base_;
    bool owns_sparse_base_ = false;
    Sets ordered_;
    std::vector<std::vector<std::size_t>> position_, outgoing_;
    std::vector<Block> blocks_;
    bool ready_ = false;
    std::uint64_t logical_bytes_ = 0;
    double prepare_ms_ = 0;
    DomainResult condition_sparse(const Required& required) const;
};
} // namespace mcviz
