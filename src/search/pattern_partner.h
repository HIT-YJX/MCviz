#pragma once
#include "pattern_core.h"
#include <cstdint>

namespace mcviz {
// Immutable adjacency over a safe global base. Each thread owns reduction masks,
// core workspace and statistics. A refused representation disables only this rule.
struct PartnerStats {
    std::uint64_t work = 0, rounds = 0, removed = 0, core_removed = 0, budget_stops = 0;
};
class PartnerRows {
public:
    PartnerRows(const Graph& g, const LayerPattern& p, const Sets& base, std::uint64_t memory);
    bool ready() const { return ready_; }
    std::uint64_t bytes() const { return bytes_; }
    CoreResult reduce(const Thresholds& beta, Sets domains, const Required& fixed,
                      std::uint64_t budget, PartnerStats& stats, CoreWorkspace& workspace) const;
private:
    struct Block { Layer target; std::size_t words; std::vector<std::uint64_t> bits; };
    const Graph& g_;
    const LayerPattern& p_;
    Sets base_;
    std::vector<std::vector<std::size_t>> pos_;
    std::vector<std::vector<Block>> blocks_;
    std::uint64_t bytes_ = 0;
    bool ready_ = false;
};
} // namespace mcviz
