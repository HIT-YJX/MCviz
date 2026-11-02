#include "contract.h"
#include "text.h"

namespace mcviz {
Thresholds parse_thresholds(const std::string& csv, unsigned) {
    Thresholds out;
    for (const auto& s : csv_fields(csv)) out.push_back(number(s, "threshold"));
    return out;
}

QueryRequest effective_query(const Graph&, const BuildSpec& spec, const QueryRequest& query) {
    auto out = query;
    if (out.gamma.empty()) out.gamma = spec.beta;
    return out;
}

Required local_required(const Graph& g, const QueryRequest& query) {
    if (!query.required) return {};
    const auto q = *query.required;
    return {{q.layer, g.get_pos(q.layer, q.id)}};
}
}
