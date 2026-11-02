#include "cli.h"
#include "text.h"
#include <algorithm>
#include <map>
#include <set>
#include <cstdlib>
#include <filesystem>

namespace mcviz {
namespace {

VertexRef parse_required(const std::string& s, unsigned K) {
    auto pos = s.find(':');
    return {static_cast<Layer>(number(s.substr(0, pos), "required layer", K - 1)), number(s.substr(pos + 1), "required id")};
}
}


// Read command options and apply the release defaults.
CliOptions parse_cli(int argc, const char* const* argv, const CliCapabilities& caps) {
    CliOptions out;
    if (caps.cohesion_build) out.pivot_memory_bytes = 512ULL << 20;
    out.index_dir = (std::filesystem::path(MCVIZ_RELEASE_ROOT) / "Index").string();
    out.log_dir = (std::filesystem::path(MCVIZ_RELEASE_ROOT) / "log").string();
    if (argc < 2 || std::string(argv[1]) == "--help") { out.help = true; return out; }
    out.command = argv[1];
    std::map<std::string, std::string> a;

    for (int i = 2; i < argc; ++i) {
        std::string key = argv[i];

        std::string val = argv[++i];
        a.emplace(key, val);
    }

    out.layers = static_cast<unsigned>(number(a.at("--layers"), "layers", 4096));

    out.graph_file = a.at("--graph");
    if (a.count("--vertices")) out.vertices_file = a.at("--vertices");

    if (a.count("--index")) out.index_file = a.at("--index");
    if (a.count("--index-dir")) out.index_dir = a.at("--index-dir");
    if (a.count("--log-dir")) out.log_dir = a.at("--log-dir");

    if (a.count("--pattern")) out.pattern_selector = a.at("--pattern");
    if (a.count("--path-order-file")) out.path_order_file = a.at("--path-order-file");
    if (a.count("--primary")) out.primary = static_cast<Layer>(number(a.at("--primary"), "primary (0-based)", out.layers - 1));
    if (a.count("--objective")) out.objective = parse_objective(a.at("--objective"));
    if (a.count("--beta")) out.beta = parse_thresholds(a.at("--beta"), out.layers);
    if (a.count("--gamma")) out.query.gamma = parse_thresholds(a.at("--gamma"), out.layers);
    if (a.count("--required")) out.query.required = parse_required(a.at("--required"), out.layers);

    if (a.count("--pivot-memory-mb")) {
        out.pivot_memory_bytes = number(a.at("--pivot-memory-mb"), "pivot memory MiB",
            std::numeric_limits<std::size_t>::max() / (1ULL << 20)) * (1ULL << 20);
    }
    if (a.count("--query-strategy")) require(a.at("--query-strategy") == "seed-core", "The release uses seed-core queries");
    if (a.count("--query-heuristic-ms")) out.query_heuristic_ms = number(a.at("--query-heuristic-ms"), "query heuristic milliseconds", 365ULL * 86400000);
    out.build_threads_explicit = a.count("--build-threads");
    if (out.build_threads_explicit) {
        out.build_threads = static_cast<unsigned>(number(a.at("--build-threads"), "build threads", 64));
    }
    if (a.count("--cohesion-work-budget"))
        out.cohesion_work_budget = number(a.at("--cohesion-work-budget"), "cohesion work budget");
    if (a.count("--cohesion-cache-mb"))
        out.cohesion_cache_bytes = number(a.at("--cohesion-cache-mb"), "cohesion cache MiB") * (1ULL << 20);
    if (a.count("--cohesion-strategy")) out.cohesion_strategy = a.at("--cohesion-strategy");
    require(out.cohesion_strategy == "auto" || out.cohesion_strategy == "general-dynamic-all", "Unsupported release construction strategy");
    if (a.count("--solver")) require(a.at("--solver") == "auto", "The release uses solver=auto");
    if (a.count("--build-core")) require(a.at("--build-core") == "queue", "The release uses queue core");

    if (a.count("--edge-bound")) require(a.at("--edge-bound") == "auto", "The release uses automatic upper bounds");
    if (a.count("--build-filter")) out.build_filter = a.at("--build-filter");
    require(out.build_filter == "auto" || out.build_filter == "off", "Unsupported build filter");

    if (a.count("--method")) {
        out.method = a.at("--method");
    }
    if (a.count("--min-layers")) {
        out.min_layers = static_cast<unsigned>(number(a.at("--min-layers"), "min-layers", 4096));
    }

    for (auto file : {&out.graph_file, &out.vertices_file}) {
        if (!file->empty() && !std::filesystem::path(*file).is_absolute() && !std::filesystem::exists(*file))
            *file = (std::filesystem::path(MCVIZ_RELEASE_ROOT) / "datasets" / *file).string();
    }
    if (!out.index_file.empty() && !std::filesystem::path(out.index_file).is_absolute())
        out.index_file = (std::filesystem::path(out.index_dir) / out.index_file).string();
    return out;
}

// Build the pattern and thresholds from command options.
BuildSpec build_spec(const CliOptions& opt) {
    BuildSpec spec;

    spec.pattern = make_pattern(opt.layers, *opt.pattern_selector, opt.primary.value_or(0), opt.path_order_file);
    spec.objective = *opt.objective;
    if (opt.beta) spec.beta = *opt.beta;
    else if (opt.command == "solve") spec.beta = opt.query.gamma;
    return spec;
}

std::string cli_help() {
    return "MCviz/MCviz+: construct/query. Baselines: solve/construct/query.\n"
           "--graph G --vertices V(optional) --layers K (2..4096)\n"
           "--pattern star|complete|path|cycle|empty|FILE|file:FILE\n"
           "--path-order-file FILE: path only, one 1-based CSV permutation.\n"
           "All other layer fields are 0-based. --primary selects the star center or optional anchor.\n"
           "--objective edge|vertex|balanced --beta positive_csv --gamma positive_csv\n"
           "--required layer:original_id --index FILE\n"
           "--index-dir DIR --log-dir DIR (defaults: release/Index and release/log).\n"
           "construct reuses compatible indexes; index queries require gamma>=beta.\n"
           "Query without --index requires --pattern/--objective/--beta for lookup.\n"
           "MCviz/MCviz+: --pivot-memory-mb N (default 512; 0 uses cover fallback).\n"
           "Construction: --build-threads 1..64 (default 32 on dynamic/support builds, otherwise 1).\n"
           "MCviz uses general-dynamic-all; queries use seed-core with --query-heuristic-ms N (default 1000).\n"
           "Baselines: --method bcvizm|motif|kpc|kspc; --min-layers K for KSPC-LP.\n"
           "Input filenames are resolved under release/datasets when absent in the current directory.\n";
}

}
