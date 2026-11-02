#pragma once
#include "pattern_core.h"

namespace mcviz {
struct RulesUpperResult {
    bool possible = true, tree_used = false, core_used = false, budget_fallback = false;
    Score upper = 0;
};
RulesUpperResult reviewed_upper(const Graph& g, const LayerPattern& p, Objective obj,
                                const Thresholds& lower, const Sets& domains, const Required& required,
                                std::uint64_t& work_budget,
                                CoreBackend core_backend = CoreBackend::reference,
                                CoreWorkspace* core_workspace = nullptr);
Thresholds improving_thresholds(const LayerPattern& p, Objective obj, const Sets& domains,
                                const Thresholds& gamma, const Sets& mandatory, Score incumbent);
}
