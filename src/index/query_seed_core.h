#pragma once
#include "pattern_index.h"
namespace mcviz {
bool seed_core_query(const Graph&, const PatternIndex&, const QueryRequest&, const Required&,
    SolverOptions, SearchResult&, SolverStats&, RunInfo&);
}
