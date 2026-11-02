#pragma once
#include "pattern_index.h"
namespace mcviz {
struct BuildFilter {
    Sets domains;
    std::vector<unsigned char> keep_edges;
    bool possible = true;
};
BuildFilter prepare_build_filter(PatternIndex& idx, const Graph& g, const BuildSpec& spec, bool possible);
void apply_build_filter(PatternIndex& idx, const Graph& g, const BuildSpec& spec, const BuildFilter& filter);
}
