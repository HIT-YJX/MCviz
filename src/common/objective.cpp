#include "objective.h"
#include "text.h"
#include <algorithm>
#include <limits>
#include <stdexcept>

namespace mcviz {

Objective parse_objective(const std::string& name) {
    if (name == "edge" || name == "MESC") return Objective::edge;
    if (name == "vertex" || name == "MVSC") return Objective::vertex;
    if (name == "balanced" || name == "MBSC") return Objective::balanced;
    throw std::runtime_error("Unknown objective: " + name);
}


const char* objective_name(Objective obj) {
    switch (obj) {
    case Objective::edge: return "edge";
    case Objective::vertex: return "vertex";
    case Objective::balanced: return "balanced";
    }
    throw std::runtime_error("Invalid objective");
}


Score checked_add(Score a, Score b) {
    require(b <= std::numeric_limits<Score>::max() - a, "Score addition overflow");
    return a + b;
}


Score checked_mul(Score a, Score b) {
    require(a == 0 || b <= std::numeric_limits<Score>::max() / a, "Score multiplication overflow");
    return a * b;
}


Score score(const LayerPattern& p, Objective obj, const std::vector<Score>& sizes) {
    require(sizes.size() == p.K && p.K >= 2, "Score dimension mismatch");
    Score value = 0;
    if (obj == Objective::edge) {
        for (auto e : p.edges) value = checked_add(value, checked_mul(sizes[e.first], sizes[e.second]));
    } else if (obj == Objective::vertex) {
        for (auto n : sizes) value = checked_add(value, n);
    } else if (obj == Objective::balanced) {
        value = *std::min_element(sizes.begin(), sizes.end());
    } else throw std::runtime_error("Invalid objective");
    return value;
}
}
