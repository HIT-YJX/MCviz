#include "pattern_index.h"
#include "common/text.h"
#include <algorithm>
#include <atomic>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
namespace mcviz {
namespace {

constexpr unsigned state_schema = 2;
const std::string state_magic = "MCVIZ-LAYER-STATE " + std::to_string(state_schema);

Score take(std::istream& in, const std::string& field, Score max = std::numeric_limits<Score>::max()) {
    std::string s;
    in >> s;

    return number(s, field, max);
}


std::string quoted(std::istream& in) {
    std::string s;
    in >> std::quoted(s);
    return s;
}



// Write a vector with its element count.
template<class T> void put_vector(std::ostream& out, const std::vector<T>& a) {
    out << a.size();
    for (auto x : a) out << ' ' << x;
    out << '\n';
}

template<class T> std::vector<T> get_vector(std::istream& in, std::size_t count_max, Score value_max) {
    auto n = take(in, "vector length", count_max);
    std::vector<T> a;
    a.reserve(static_cast<std::size_t>(n));
    for (Score j = 0; j < n; ++j) a.push_back(static_cast<T>(take(in, "vector entry", value_max)));

    return a;
}
}


void PatternIndex::save(const Graph& g, const std::string& path) const {
    std::ostringstream out;

    const auto& spec = identity.spec;
    const auto& p = spec.pattern;
    out << identity.schema << ' ' << std::quoted(identity.contract) << ' ' << std::quoted("pattern-engine-v2") << '\n';
    out << std::quoted(identity.engine) << ' ' << std::quoted(identity.value_kind) << ' '
        << std::quoted(identity.graph_hash) << ' ' << std::quoted(identity.pattern_hash) << '\n';
    const auto solver = options.solver;
    out << std::quoted(query_execution_path()) << ' ' << std::quoted(bound) << ' ' << std::quoted(solver) << '\n';
    out << std::quoted(source_graph_sha256) << ' ' << std::quoted(source_vertices_sha256) << '\n';
    out << p.K << ' ' << std::quoted(objective_name(spec.objective)) << '\n';

    put_vector(out, spec.beta);

    out << p.edges.size() << '\n';
    for (auto e : p.edges) out << e.first << ' ' << e.second << '\n';

    out << std::quoted(p.selector) << ' ' << std::quoted(p.source) << ' ' << p.anchor_layer << '\n';
    put_vector(out, p.path_order_input); put_vector(out, p.path_order_internal);
    out << "vertex_ids " << g.K << '\n';

    for (const auto& ids : g.vertex_ids) put_vector(out, ids);

    for (const auto& c : core_domains) put_vector(out, c);
    out << weights.size() << '\n';
    for (const auto& w : weights) out << w.edge.u.layer << ' ' << w.edge.u.pos << ' '
        << w.edge.v.layer << ' ' << w.edge.v.pos << ' ' << w.weight << '\n';

    put_vector(out, forest.edge_ids);
    out << forest.nodes.size() << '\n';
    for (const auto& n : forest.nodes) { out << n.height << ' '; put_vector(out, n.children); }
    put_vector(out, forest.roots);
    const auto payload = out.str();
    const auto header = state_magic + "\n" + std::to_string(payload.size()) + " " + sha256(payload) + "\n";
    std::ofstream file(path, std::ios::binary);
    file << header << payload;
}


PatternIndex PatternIndex::load(const Graph& g, const std::string& path, Engine engine) {
    std::ifstream file(path, std::ios::binary | std::ios::ate);
    file.seekg(0);
    std::string line;
    std::getline(file, line);
    std::getline(file, line);
    auto header = fields(line);
    auto length = number(header[0], "payload length", 512ULL * 1024 * 1024);
    std::string payload(static_cast<std::size_t>(length), '\0');

    file.read(payload.data(), static_cast<std::streamsize>(length));

    std::istringstream in(payload);

    PatternIndex idx;
    idx.engine = engine;
    auto& id = idx.identity;
    id.schema = static_cast<unsigned>(take(in, "schema version", schema_version));
    id.contract = quoted(in);
    quoted(in);
    id.engine = quoted(in); id.value_kind = quoted(in); id.graph_hash = quoted(in); id.pattern_hash = quoted(in);
    idx.execution_path = quoted(in); idx.bound = quoted(in); idx.options.solver = quoted(in);
    idx.options.edge_bound = idx.bound;
    idx.source_graph_sha256 = quoted(in); idx.source_vertices_sha256 = quoted(in);
    auto& spec = id.spec;
    auto& p = spec.pattern;

    p.K = static_cast<unsigned>(take(in, "K", 4096));
    spec.objective = parse_objective(quoted(in));
    spec.beta = get_vector<Score>(in, g.K, std::numeric_limits<Score>::max());
    auto edges = take(in, "pattern edge count", Score(g.K) * (g.K - 1) / 2);
    for (Score j = 0; j < edges; ++j) {
        Layer a = static_cast<Layer>(take(in, "pattern layer", g.K - 1));
        Layer b = static_cast<Layer>(take(in, "pattern layer", g.K - 1));

        p.edges.emplace_back(a, b);
    }
    p.selector = quoted(in); p.source = quoted(in); p.anchor_layer = static_cast<Layer>(take(in, "anchor", g.K - 1));
    p.path_order_input = get_vector<Layer>(in, g.K, g.K);
    p.path_order_internal = get_vector<Layer>(in, g.K, g.K - 1);
    p.normalize();
    std::string domain_tag;
    in >> domain_tag;
    take(in, "vertex ID layer count", g.K);

    for (Layer i = 0; i < g.K; ++i) {
        auto ids = get_vector<VertexId>(in, g.vertex_ids[i].size(), std::numeric_limits<VertexId>::max());
    }

    for (Layer i = 0; i < g.K; ++i) {
        auto c = get_vector<VertexPos>(in, g.vertex_ids[i].size(), g.vertex_ids[i].empty() ? 0 : g.vertex_ids[i].size() - 1);
        idx.core_domains.push_back(std::move(c));
    }
    std::vector<GraphEdge> relevant;
    if (p.connected()) for (auto e : g.edges) if (p.adjacent(e.u.layer, e.v.layer)) relevant.push_back(e);

    auto nw = take(in, "weight count", relevant.size());
    for (auto expected : relevant) {
        take(in, "edge layer"); take(in, "edge pos"); take(in, "edge layer"); take(in, "edge pos");

        idx.weights.push_back({expected, take(in, "edge weight")});
    }
    auto& f = idx.forest;
    if (p.connected()) {
        f.offsets.push_back(0);

        for (Layer i = 0; i < g.K; ++i) {
            for (VertexPos u = 0; u < g.vertex_ids[i].size(); ++u) f.vertices.push_back({i, u});
            f.offsets.push_back(f.vertices.size());
        }
    }
    std::size_t n = f.vertices.size();
    f.edge_ids = get_vector<std::size_t>(in, n, nw == 0 ? 0 : nw - 1);
    auto nn = take(in, "hierarchy node count", 2 * n);
    for (Score j = 0; j < nn; ++j) {
        Score h = take(in, "merge height");
        auto children = get_vector<std::size_t>(in, n, nn == 0 ? 0 : nn - 1);
        f.nodes.push_back({h, std::move(children)});
    }

    f.roots = get_vector<std::size_t>(in, n, nn == 0 ? 0 : nn - 1);
    in >> std::ws;
    return idx;
}
}
