#include "api.h"
#include "bcvizm_lp.h"
#include "../common/cli.h"
#include "../common/index_library.h"
#include "../common/log.h"
#include <chrono>
#include <filesystem>
#include <iostream>

namespace {
using Clock = std::chrono::steady_clock;

double seconds(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}
}

int main(int argc, char** argv) {
    const auto began = Clock::now();
    mcviz::CliCapabilities caps;

    caps.commands = {"solve", "construct", "query"};
    caps.method = true;
    caps.min_layers = true;
    auto opt = mcviz::parse_cli(argc, argv, caps);
    if (opt.help) {
std::cout << mcviz::cli_help();
        return 0;
    }
    const auto g = mcviz::Graph::load(opt.layers, opt.graph_file, opt.vertices_file);

    mcviz::BuildSpec spec;
    mcviz::baseline::BcvizmIndex index;
    auto path = opt.index_file;

    bool reused = false;

    double build_s = 0;

    if (opt.command == "query") {
        if (path.empty()) {
            spec = mcviz::build_spec(opt);
            const auto paths = mcviz::library_paths(opt, g, spec, "BCviz-M");
            path = mcviz::find_library_index(paths, g, spec, "BCviz-M-LP", "safe-upper");
        }
        index = mcviz::baseline::load_bcvizm(path);
        spec = index.identity.spec;
    } else {
        spec = mcviz::build_spec(opt);
        if (opt.method == "bcvizm") {
            if (opt.command == "construct") {
                const auto paths = mcviz::library_paths(opt, g, spec, "BCviz-M");

                auto existing = !path.empty() && std::filesystem::exists(path) ? path
                    : mcviz::find_library_index(paths, g, spec, "BCviz-M-LP", "safe-upper");
                if (!existing.empty()) {
                    index = mcviz::baseline::load_bcvizm(existing);
                    path = existing;
                    reused = true;
                } else {
                    path = opt.index_file.empty() ? paths.index : opt.index_file;
                }
            }
            if (!reused) {
                const auto start = Clock::now();
                index = mcviz::baseline::construct_bcvizm(g, spec);
                build_s = seconds(start);
            }
        }
    }

    if (opt.command == "construct") {
        if (!reused) {
            std::filesystem::create_directories(std::filesystem::path(path).parent_path());
            mcviz::baseline::save_bcvizm(index, path);
        }
        if (reused && !opt.index_file.empty() && opt.index_file != path) {
            std::filesystem::create_directories(std::filesystem::path(opt.index_file).parent_path());
            std::filesystem::copy_file(path, opt.index_file);
            path = opt.index_file;
        }
        mcviz::log_build(opt, spec, opt.method, build_s, std::filesystem::file_size(path), reused);

        return 0;
    }

    const auto query = mcviz::effective_query(g, spec, opt.query);
    const auto result = opt.method == "bcvizm"
        ? mcviz::baseline::query_bcvizm(g, index, query)
        : mcviz::baseline::run(g, spec.pattern, spec.objective,
                               query, opt.method, opt.min_layers.value_or(0));
    const auto elapsed = seconds(began);
    mcviz::log_search(opt, spec, query, result, opt.method, elapsed);
    return 0;
}
