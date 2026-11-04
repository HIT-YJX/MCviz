#include "engine_cli.h"
#include "common/cli.h"
#include "common/index_library.h"
#include "common/log.h"
#include "common/text.h"
#include <chrono>
#include <filesystem>
#include <iostream>

namespace mcviz {
namespace {
using Clock = std::chrono::steady_clock;
double seconds(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}
}

int run_engine(int argc, char** argv, Engine engine) {
    const auto begin = Clock::now();

    CliCapabilities caps;
    caps.commands = {"construct", "query"};
    caps.cohesion_build = caps.solver = true;
    caps.support_build = engine == Engine::upper;
    auto opt = parse_cli(argc, argv, caps);
    if (opt.help) {
std::cout << cli_help();
        return 0;
    }
    require(engine == Engine::exact || opt.cohesion_strategy == "auto", "MCviz+ uses automatic upper-bound construction");
    const auto g = Graph::load(opt.layers, opt.graph_file, opt.vertices_file);
    const std::string method = engine == Engine::exact ? "mcviz" : "mcviz-plus";
    const std::string name = engine == Engine::exact ? "MCviz" : "MCviz-plus";
    PatternIndex idx;
    auto path = opt.index_file;

    double build_seconds = 0;
    bool reused = false;
    if (opt.command == "construct") {
        const auto spec = build_spec(opt);
        require(opt.cohesion_strategy == "auto" ||
                (engine == Engine::exact && spec.objective == Objective::edge &&
                 !spec.pattern.edges.empty() && spec.pattern.connected()),
                "general-dynamic-all requires a connected exact edge build");
        if (!opt.build_threads_explicit) {
            const auto& p = spec.pattern;
            bool dynamic = engine == Engine::exact && spec.objective == Objective::edge && !p.edges.empty() && p.connected();
            bool support = engine == Engine::upper && !p.star_centers().empty();

            if (dynamic || support) opt.build_threads = 32;
        }
        const auto paths = library_paths(opt, g, spec, name);
        auto existing = !path.empty() && std::filesystem::exists(path) ? path
            : find_library_index(paths, g, spec, engine_name(engine));
        if (!existing.empty()) {
            idx = PatternIndex::load(g, existing, engine);
            path = existing;
            reused = true;
        } else {
            IndexOptions options;
            options.search.pivot_memory_bytes = opt.pivot_memory_bytes;
            options.build_threads = opt.build_threads;

            options.build_filter = opt.build_filter;
            options.cohesion_work_budget = opt.cohesion_work_budget;
            options.cohesion_cache_bytes = opt.cohesion_cache_bytes;
            const auto start = Clock::now();
            idx = PatternIndex::construct(g, spec, engine, options);

            build_seconds = seconds(start);
            path = opt.index_file.empty() ? paths.index : opt.index_file;

            std::filesystem::create_directories(std::filesystem::path(path).parent_path());
            idx.save(g, path);
        }
        if (reused && !opt.index_file.empty() && opt.index_file != path) {
            std::filesystem::create_directories(std::filesystem::path(opt.index_file).parent_path());
            std::filesystem::copy_file(path, opt.index_file);
            path = opt.index_file;
        }
        if (!reused && (idx.build_counters["build_shared_budget_fallback"] ||
                        idx.build_counters["build_shared_memory_fallback"] ||
                        idx.build_counters["build_anchor_budget_fallback"] ||
                        idx.build_counters["build_anchor_memory_fallback"])) {
std::cout << "star construction status=exact-fallback\n";
        }
        log_build(opt, idx.identity.spec, method, build_seconds, std::filesystem::file_size(path), reused);
        return 0;
    }

    if (path.empty()) {
        const auto spec = build_spec(opt);

        const auto paths = library_paths(opt, g, spec, name);
        path = find_library_index(paths, g, spec, engine_name(engine));
    }
    idx = PatternIndex::load(g, path, engine);

    idx.options.search.pivot_memory_bytes = opt.pivot_memory_bytes;
    idx.options.query_heuristic_ms = opt.query_heuristic_ms;
    RunInfo run;

    const auto q = effective_query(g, idx.identity.spec, opt.query);
    const auto result = idx.query(g, q, run);
    const auto elapsed = seconds(begin);
    log_search(opt, idx.identity.spec, q, result, method, elapsed);
    return result.status == Status::optimal || result.status == Status::infeasible ? 0 : 1;
}
}
