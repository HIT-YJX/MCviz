#pragma once
#include "pattern_reuse_build.h"
namespace mcviz {
// Exact edge weights with dynamic scheduling, shared completed-domain scores,
// partner reduction and witnessed lower bounds. Search remains exponential.
PatternReuseBuildResult build_pattern_dynamic(const Graph& g, const LayerPattern& pattern,
    const Thresholds& beta, Sets& base, unsigned threads, std::uint64_t cache_bytes,
    std::uint64_t partner_work, std::uint64_t pivot_memory_bytes = 64ULL << 20);
}
