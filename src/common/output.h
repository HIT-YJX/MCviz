#pragma once
#include "contract.h"
#include <map>

namespace mcviz {
struct RunInfo {
    std::string algorithm, execution_path, solver, bound;
    bool accelerated = false;
    std::uint64_t regions = 0, solver_calls = 0;
    std::map<std::string, double> timing;
    std::map<std::string, std::uint64_t> counters;
};

}
