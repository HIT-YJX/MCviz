#include "bcvizm_lp.h"
#include "../common/text.h"
#include <fstream>
#include <iomanip>
#include <sstream>

namespace mcviz::baseline {
namespace {
// Serialize index identity, typed vertices, and the search order.
std::string payload(const BcvizmIndex& index) {
    const auto& id = index.identity;
    const auto& P = id.spec.pattern;
    std::ostringstream out;
    out << "BCVIZM_LP_2 " << id.schema << ' ' << id.contract << ' ' << P.K << '\n';
    out << id.engine << ' ' << id.value_kind << ' ' << id.graph_hash << ' '
        << id.pattern_hash << '\n';
    out << objective_name(id.spec.objective) << ' ' << P.anchor_layer << '\n';

    out << std::quoted(P.selector) << ' ' << std::quoted(P.source) << '\n';

    out << P.path_order_input.size();

    for (auto p : P.path_order_input) out << ' ' << p;
    out << '\n' << P.path_order_internal.size();
    for (auto p : P.path_order_internal) out << ' ' << p;
    out << '\n' << id.spec.beta.size();
    for (auto t : id.spec.beta) out << ' ' << t;

    out << '\n' << P.edges.size() << '\n';
    for (auto [i,j] : P.edges) out << i << ' ' << j << '\n';
    out << index.domain.size() << '\n';
    for (auto v : index.domain) out << v.layer << ' ' << v.id << '\n';
    out << index.order.size() << '\n';
    for (auto e : index.order) out << e.vertex << ' ' << e.key << ' ' << e.root << '\n';
    return out.str();
}

}

void save_bcvizm(const BcvizmIndex& index, const std::string& path) {
    const auto body = payload(index);
    std::ofstream out(path, std::ios::binary);
    out << body << "CHECKSUM " << sha256(body) << '\n';
}

BcvizmIndex load_bcvizm(const std::string& path) {
    const auto bytes = read_bytes(path);
    const auto mark = bytes.rfind("CHECKSUM ");
    std::istringstream in(bytes.substr(0, mark));

    BcvizmIndex index;
    auto& id = index.identity;
    auto& P = id.spec.pattern;
    std::string magic, objective;
    in >> magic >> id.schema >> id.contract >> P.K;
    in >> id.engine >> id.value_kind >> id.graph_hash >> id.pattern_hash;

    in >> objective >> P.anchor_layer >> std::quoted(P.selector) >> std::quoted(P.source);
    id.spec.objective = parse_objective(objective);

    size_t count = 0;
    in >> count;
    P.path_order_input.resize(count);

    for (auto& p : P.path_order_input) in >> p;
    in >> count;
    P.path_order_internal.resize(count);

    for (auto& p : P.path_order_internal) in >> p;
    in >> count;
    id.spec.beta.resize(count);
    for (auto& t : id.spec.beta) in >> t;
    in >> count;

    P.edges.resize(count);
    for (auto& e : P.edges) in >> e.first >> e.second;
    P.normalize();
    in >> count;
    index.domain.resize(count);
    for (auto& v : index.domain) in >> v.layer >> v.id;
    in >> count;
    index.order.resize(count);
    for (auto& e : index.order) {
        unsigned root = 0;
        in >> e.vertex >> e.key >> root;
        e.root = root != 0;
    }
    return index;
}


// Load the saved BCviz-M index through the shared file reader.
BcvizmIndex load_bcvizm(const Graph&, const std::string& path) {
    return load_bcvizm(path);
}
}
