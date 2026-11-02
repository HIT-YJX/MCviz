#include "index_library.h"
#include "text.h"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>

namespace mcviz {
namespace {
namespace fs = std::filesystem;
// Map dataset filenames to their index naming aliases.
std::string alias(const std::string& stem) {
    static const std::map<std::string, std::string> names{
        {"imdb_hgb_k4","imdb"},{"dblp_hgb_k4","dblp"},{"DBLP_k2","dblp2"},
        {"lastfm_hetrec_k3","lastfm"},{"acm_k4","acm"},{"bibsonomy_k3","bibsonomy"},
        {"aminer_k3","aminer"},{"WikiPedia_k2","wikipedia"},{"marvel_k2","marvel"},
        {"movielens25m_k3","movielens"},{"writer_k2","writer"},{"youtube","youtube"},
        {"hetionet_k6","hetionet"},{"drkg_k5","drkg"},{"wordnet_k25","wordnet"}};
    auto i = names.find(stem);
    auto name = i == names.end() ? stem : i->second;
    return name;
}



std::string token(std::istream& in) { std::string s; in >> s; return s; }

std::string quoted(std::istream& in) { std::string s; in >> std::quoted(s); return s; }


// Read a bounded integer from index metadata.
Score take(std::istream& in, Score limit) { return number(token(in), "index metadata", limit); }

// Read a length-prefixed threshold vector.
Thresholds vector(std::istream& in, std::size_t limit) {
    auto n = take(in, limit); Thresholds a;

    for (Score j = 0; j < n; ++j) a.push_back(take(in, std::numeric_limits<Score>::max()));

    return a;
}


// Read pattern edges and normalize their order.
void edges(std::istream& in, LayerPattern& p) {
    auto n = take(in, Score(p.K) * (p.K - 1) / 2);
    for (Score j = 0; j < n; ++j) {
        auto u = static_cast<Layer>(take(in, p.K - 1));
        auto v = static_cast<Layer>(take(in, p.K - 1));
        p.edges.emplace_back(u, v);
    }
    p.normalize();
}


// Recognize existing paths, including symbolic links.
bool exists(const std::string& path) { return fs::exists(path) || fs::is_symlink(path); }
}


// Derive index paths from the dataset and query configuration.
LibraryPaths library_paths(const CliOptions& opt, const Graph& g, const BuildSpec& spec,
                           const std::string& algorithm) {
    const auto& p = spec.pattern;

    auto shape = p.selector;
    std::string pattern_suffix;
    if (shape != "star" && shape != "complete" && shape != "path" && shape != "cycle" && shape != "empty")
        { shape = "custom"; pattern_suffix = "-ph" + p.hash().substr(0, 12); }
    if (shape == "path") {
        auto normal = make_pattern(g.K, "path");
        if (normal.edges != p.edges) pattern_suffix = "-ph" + p.hash().substr(0, 12);
    }

    auto b = *std::max_element(spec.beta.begin(), spec.beta.end());

    bool standard = true, all = true, singleton = false;
    for (Layer j = 0; j < g.K; ++j) {
        singleton = singleton || g.vertex_ids[j].size() == 1;
        all = all && spec.beta[j] == b;
        standard = standard && spec.beta[j] == (g.vertex_ids[j].size() == 1 ? 1 : b);
    }
    bool simple = all || standard;

    std::string beta = "b" + std::to_string(b);

    if (!simple) {
        beta = "bvec";
        for (Layer j = 0; j < g.K; ++j) { if (j) beta += 'x'; beta += std::to_string(spec.beta[j]); }
    }

    auto prefix = alias(fs::path(opt.graph_file).stem().string()) + '_' + algorithm + '_';
    auto stem = prefix + shape + '_' + objective_name(spec.objective) + '_' + beta;

    if (shape == "star") stem += "_p" + std::to_string(p.anchor_layer);
    stem += pattern_suffix;
    auto base = (fs::path(opt.index_dir) / (stem + ".idx")).string();
    if (simple && singleton && b > 1 && exists(base)) {
        auto prior = library_index_identity(base);
        if (prior.graph_hash == g.hash() && prior.pattern_hash == p.hash() &&
            prior.spec.objective == spec.objective && prior.spec.beta != spec.beta) {
            if (all) stem += "-all";
            else {
                stem += "-bv";
                for (Layer j = 0; j < g.K; ++j) { if (j) stem += 'x'; stem += std::to_string(spec.beta[j]); }
            }
        }
    }
    return {(fs::path(opt.index_dir) / (stem + ".idx")).string(),
            stem, prefix};
}



BuildIdentity library_index_identity(const std::string& path) {
    std::ifstream in(path, std::ios::binary);

    auto magic = token(in);
    BuildIdentity id;
    auto& p = id.spec.pattern;
    if (magic == "MCVIZ-LAYER-STATE") {
        take(in, 2);
        take(in, 512ULL * 1024 * 1024); token(in);
        id.schema = static_cast<unsigned>(take(in, schema_version)); id.contract = quoted(in);
        quoted(in);
        id.engine = quoted(in); id.value_kind = quoted(in); id.graph_hash = quoted(in); id.pattern_hash = quoted(in);
        quoted(in); quoted(in); quoted(in); quoted(in); quoted(in);
        p.K = static_cast<unsigned>(take(in, 4096));
        id.spec.objective = parse_objective(quoted(in)); id.spec.beta = vector(in, p.K);
        edges(in, p); p.selector = quoted(in); p.source = quoted(in); p.anchor_layer = static_cast<Layer>(take(in, p.K - 1));
    } else if (magic == "BCVIZM_LP_2") {
        id.schema = static_cast<unsigned>(take(in, schema_version)); id.contract = token(in);
        p.K = static_cast<unsigned>(take(in, 4096));
        id.engine = token(in); id.value_kind = token(in); id.graph_hash = token(in); id.pattern_hash = token(in);

        id.spec.objective = parse_objective(token(in)); p.anchor_layer = static_cast<Layer>(take(in, p.K - 1));

        p.selector = quoted(in); p.source = quoted(in);
        vector(in, p.K); vector(in, p.K); id.spec.beta = vector(in, p.K); edges(in, p);
    } else throw std::runtime_error("Unsupported public index metadata format: " + path);
    return id;
}

std::string find_library_index(const LibraryPaths& paths, const Graph& g, const BuildSpec& spec,
                               const std::string& engine, const std::string& value_kind) {
    auto dir = fs::path(paths.index).parent_path();
    if (!fs::exists(dir)) return {};

    std::vector<std::string> candidates;
    if (exists(paths.index)) candidates.push_back(paths.index);
    for (const auto& e : fs::directory_iterator(dir)) {
        auto name = e.path().filename().string();
        if (e.path().extension() == ".idx" && name.rfind(paths.prefix, 0) == 0 && e.path().string() != paths.index)
            candidates.push_back(e.path().string());
    }
    const auto graph_hash = g.hash(), pattern_hash = spec.pattern.hash();
    for (const auto& path : candidates) {
        auto id = library_index_identity(path);
        if (id.graph_hash == graph_hash && id.pattern_hash == pattern_hash && id.engine == engine &&
            id.spec.objective == spec.objective && id.spec.beta == spec.beta &&
            (spec.pattern.selector != "star" || id.spec.pattern.anchor_layer == spec.pattern.anchor_layer) &&
            (value_kind.empty() || id.value_kind == value_kind)) return path;
    }

    return {};
}

}
