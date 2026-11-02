#pragma once
#include "common/contract.h"

namespace mcviz {
struct DomainResult { Sets domains; bool possible = true; };
DomainResult condition_domains(const Graph& g, const LayerPattern& p, Sets domains,
                               const Required& required);
bool contains_required(const Sets& domains, const Required& required);
bool enough_capacity(const Sets& domains, const Thresholds& gamma);
Score capacity_bound(const LayerPattern& p, Objective obj, const Sets& domains);
}
