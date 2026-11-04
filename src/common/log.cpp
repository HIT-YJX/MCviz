#include "log.h"
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <locale>
#include <map>
#include <sstream>

namespace mcviz {
namespace {
struct Cell {
    const char* name;
    const char* type;
    std::string value;
};

std::string lower(std::string s) {
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}



std::string dataset(const CliOptions& opt, const BuildSpec& spec) {
    auto name = lower(std::filesystem::path(opt.graph_file).stem().string());

    static const std::map<std::string, std::string> names = {
        {"wordnet_k25", "wordnet"}, {"wordnet_p3", "wordnet-p3"},
        {"lastfm_hetrec_k3", "lastfm"}, {"last.fm", "lastfm"},
        {"hetionet_k6", "hetionet"}, {"drkg_k5", "drkg"},
        {"freebase_k6", "freebase"}, {"bibsonomy_k3", "bibsonomy"},
        {"aminer_k3", "aminer"}, {"aminer_ap", "aminer-ap"},
        {"youtube_k32", "youtube"}, {"dblp_hgb_k4", "dblp"}
    };
    auto it = names.find(name);
    if (it != names.end()) name = it->second;
    const std::vector<std::pair<Layer, Layer>> tree{{0, 2}, {0, 3}, {0, 4}, {1, 4}};
    if (name == "drkg" && spec.pattern.K == 5 && spec.pattern.edges == tree) name = "drkg-tree";
    if (name == "wordnet") name = "WordNet";
    if (name == "wordnet-p3") name = "WordNet-P3";
    if (name == "lastfm") name = "LastFM";
    return name;
}

// Normalize algorithm names for timing records.
std::string method_key(const std::string& method) {
    auto name = lower(method);
    static const std::map<std::string, std::string> names = {
        {"mcviz-lp", "mcviz"}, {"mcviz+", "mcviz-plus"}, {"mcviz+-lp", "mcviz-plus"},
        {"bcviz-m", "bcvizm"}, {"bcviz-m-lp", "bcvizm"}, {"bcvizm-lp", "bcvizm"},
        {"motif-lp", "motif"}, {"kpc-lp", "kpc"}, {"kspc-lp", "kspc"},
        {"kspc-lp-all", "kspc"}
    };
    auto it = names.find(name);
    return it == names.end() ? name : it->second;
}


std::string text(const std::string& s) {
    return !s.empty() && s.front() == '\\' ? "\\" + s : s;
}

template<class T> std::string number(T value) {
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out << std::setprecision(17) << value;
    return out.str();
}


// Format thresholds as a JSON array.
std::string array(const Thresholds& a) {
    std::string s = "[";
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (i) s += ',';
        s += number(a[i]);
    }
    return s + ']';
}


std::string pattern(const LayerPattern& p) {
    std::string s = "[";
    for (std::size_t i = 0; i < p.edges.size(); ++i) {
        if (i) s += ',';
        s += '[' + number(p.edges[i].first) + ',' + number(p.edges[i].second) + ']';
    }
    return s + ']';
}



void cell(std::ostream& out, const std::string& s) {
    if (s.find_first_of("\t\r\n\"") == std::string::npos) { out << s; return; }
    out << '"';
    for (char c : s) {
        if (c == '"') out << '"';
        out << c;
    }
    out << '"';
}



// Append one typed TSV row; O(row bytes) time and space, with one schema per file.
void append(const CliOptions& opt, const char* name, const std::vector<Cell>& row) {
    std::filesystem::create_directories(opt.log_dir);
    std::ofstream out(std::filesystem::path(opt.log_dir) / name, std::ios::app);
    if (out.tellp() == 0) {
        out << "# MCviz-test2 TSV v1\n";
        for (std::size_t i = 0; i < row.size(); ++i) {
            if (i) out << '\t';

            out << row[i].name;
        }
        out << '\n';
        for (std::size_t i = 0; i < row.size(); ++i) {
            if (i) out << '\t';

            out << row[i].type;
        }
        out << '\n';
    }
    for (std::size_t i = 0; i < row.size(); ++i) {
        if (i) out << '\t';

        cell(out, row[i].value);
    }
    out << '\n';
}
}

// Record index construction time, size, and reuse status.
void log_build(const CliOptions& opt, const BuildSpec& spec, const std::string& method,
               double seconds, std::uintmax_t bytes, bool reused) {
    const auto key = method_key(method);

    const std::string elapsed = reused ? "\\N" : number(seconds);
    append(opt, "builds.txt", {
        {"dataset", "text", text(dataset(opt, spec))}, {"method", "text", text(key)},
        {"status", "text", reused ? "reused" : "ok"}, {"seconds", "number", elapsed},
        {"bytes", "number", number(bytes)}, {"value_seconds", "number", elapsed},
        {"value_mib", "number", number(static_cast<double>(bytes) / 1048576.0)},
        {"objective", "text", objective_name(spec.objective)}, {"beta", "json", array(spec.beta)},
        {"pattern", "json", pattern(spec.pattern)}, {"primary", "number", number(spec.pattern.anchor_layer)},
        {"threads", "number", number(opt.build_threads)}, {"reused", "bool", reused ? "true" : "false"}
    });
std::cout << key << " construct status=" << (reused ? "reused" : "ok") << " bytes=" << bytes;
    if (!reused) {
std::cout << " seconds=" << number(seconds);
    }
std::cout << '\n';
}


void log_search(const CliOptions& opt, const BuildSpec& spec, const QueryRequest& query,
                const SearchResult& result, const std::string& method, double seconds) {
    const auto key = method_key(method);
    const bool complete = result.status == Status::optimal || result.status == Status::infeasible;
    const auto& gamma = query.gamma.empty() ? spec.beta : query.gamma;
    const std::string required = query.required
        ? number(query.required->layer) + ':' + number(query.required->id) : "\\N";
    append(opt, "searches.txt", {
        {"dataset", "text", text(dataset(opt, spec))}, {"method", "text", text(key)},
        {"status", "text", complete ? "ok" : status_name(result.status)},
        {"found", "bool", result.found ? "true" : "false"},
        {"score", "number", result.found ? number(result.score) : "\\N"},
        {"seconds", "number", number(seconds)}, {"objective", "text", objective_name(spec.objective)},
        {"beta", "json", array(spec.beta)}, {"gamma", "json", array(gamma)},
        {"required", "text", required}, {"pattern", "json", pattern(spec.pattern)},
        {"primary", "number", number(spec.pattern.anchor_layer)}
    });

std::cout << key << " search status=" << status_name(result.status)
    << " found=" << (result.found ? "true" : "false");
    if (result.found) {
std::cout << " score=" << result.score;
    }
std::cout << " seconds=" << number(seconds) << '\n';
}



void log_conversion(const CliOptions& opt, const BuildSpec& spec, const std::string& method,
                    double forward, double reverse, std::uint64_t tree_bytes,
                    std::uint64_t linear_bytes) {
    const auto key = method_key(method);
    for (const auto& phase : {std::make_pair("tree_to_linear", forward),
                              std::make_pair("linear_to_tree", reverse)}) {
        append(opt, "conversions.txt", {
            {"dataset", "text", text(dataset(opt, spec))}, {"method", "text", text(key)},
            {"status", "text", "ok"}, {"direction", "text", phase.first},
            {"seconds", "number", number(phase.second)}, {"value_ms", "number", number(phase.second * 1000.0)},
            {"tree_bytes", "number", number(tree_bytes)}, {"linear_bytes", "number", number(linear_bytes)},
            {"objective", "text", objective_name(spec.objective)}, {"beta", "json", array(spec.beta)},
            {"pattern", "json", pattern(spec.pattern)}, {"primary", "number", number(spec.pattern.anchor_layer)}
        });
    }
std::cout << key << " conversion status=ok tree_to_linear_seconds=" << number(forward)
    << " linear_to_tree_seconds=" << number(reverse) << '\n';
}
}
