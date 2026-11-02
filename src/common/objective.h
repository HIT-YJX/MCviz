#pragma once
#include "pattern.h"

namespace mcviz {
enum class Objective { edge, vertex, balanced };
Objective parse_objective(const std::string& name);
const char* objective_name(Objective objective);
Score checked_add(Score a, Score b);
Score checked_mul(Score a, Score b);
Score score(const LayerPattern& pattern, Objective objective, const std::vector<Score>& sizes);
}
