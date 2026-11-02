#pragma once
#include "common/contract.h"

namespace mcviz {
struct PatternFilterResult {
    std::vector<unsigned char> keep_edges;
    bool mask_available = false;
    std::uint64_t eligible_pairs = 0, completed_pairs = 0;
    std::uint64_t skipped_work = 0, skipped_memory = 0, skipped_threshold = 0;
    std::uint64_t removed_edges = 0, work = 0, wedges = 0, peak_bytes = 0;
    double prepare_ms = 0, filter_ms = 0;
};

PatternFilterResult filter_butterflies(const Graph& g, const LayerPattern& p,
                                      const Thresholds& beta, const Sets& core,
                                      std::uint64_t work_budget, std::uint64_t memory_budget);
}
