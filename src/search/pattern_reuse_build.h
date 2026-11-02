#pragma once
#include "pattern_solver.h"
#include <map>
namespace mcviz {
struct PatternReuseBuildResult {
    std::vector<Score> values; // Full g.edges order; non-pattern edges remain zero.
    std::string value_kind;
    std::map<std::string, std::uint64_t> counters;
    std::map<std::string, double> timing;
    // Accounted representation/cache/array payload, excluding original Graph,
    // counter maps, allocator metadata and exact solver
    // workspace. This is not a process RSS limit or a complete peak-RSS estimate.
    std::uint64_t logical_payload_bytes = 0;
};
}
