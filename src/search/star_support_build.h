#pragma once
#include "star_solver.h"
#include <map>

namespace mcviz {
struct SupportBuildResult {
    std::vector<Score> values;
    std::string strategy;
    std::map<std::string, std::uint64_t> counters;
    std::map<std::string, double> timing;
};

SupportBuildResult build_star_support(const Graph& active, const LayerPattern& pattern,
                                     Objective objective, const Thresholds& beta,
                                     unsigned threads, std::uint64_t cache_bytes);
}
