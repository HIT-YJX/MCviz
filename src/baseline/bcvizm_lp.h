#pragma once

#include "api.h"
#include "../common/contract.h"

namespace mcviz::baseline {
struct BcvizmEntry { unsigned vertex = 0; Score key = 0; bool root = false; };
struct BcvizmIndex {
    BuildIdentity identity;
    std::vector<VertexRef> domain;
    std::vector<BcvizmEntry> order;
};
struct BcvizmStats {
    std::uint64_t regions = 0, zero_cuts = 0;
    std::uint64_t capacity_skipped = 0, stack_states = 0;
    std::uint64_t star_regions = 0, generic_regions = 0;
    std::uint64_t star_states = 0, generic_states = 0;
    std::uint64_t star_branch_states = 0, generic_branch_states = 0;
    std::uint64_t star_max_center_vertices = 0, star_max_secondary_vertices = 0;
    std::uint64_t star_saturated_leaves = 0, star_common_pruned = 0;
    std::uint64_t active_region_vertices = 0;
};

Score safe_edge_cap(const Graph& graph, const LayerPattern& pattern,
                    Objective objective, LocalVertex u, LocalVertex v);
BcvizmIndex construct_bcvizm(const Graph& graph, const BuildSpec& spec);
SearchResult query_bcvizm(const Graph& graph, const BcvizmIndex& index,
                          const QueryRequest& request, BcvizmStats* stats = nullptr);
void save_bcvizm(const BcvizmIndex& index, const std::string& path);
BcvizmIndex load_bcvizm(const std::string& path);
BcvizmIndex load_bcvizm(const Graph& graph, const std::string& path);
SearchResult run_bcvizm(const Graph& graph, const LayerPattern& pattern,
                        Objective objective, const QueryRequest& request,
                        BcvizmStats* stats = nullptr);
}
