#pragma once
#include "graph.h"
#include "objective.h"

namespace mcviz {
inline constexpr unsigned schema_version = 1;
inline constexpr const char* contract_version = "all-layer-positive-v1";

struct BuildSpec {
    LayerPattern pattern;
    Objective objective = Objective::edge;
    Thresholds beta;
};

struct BuildIdentity {
    unsigned schema = schema_version;
    std::string contract = contract_version;
    std::string engine, value_kind, graph_hash, pattern_hash;
    BuildSpec spec;
};

Thresholds parse_thresholds(const std::string& csv, unsigned K);
QueryRequest effective_query(const Graph& g, const BuildSpec& spec, const QueryRequest& query);
Required local_required(const Graph& g, const QueryRequest& query);
}
