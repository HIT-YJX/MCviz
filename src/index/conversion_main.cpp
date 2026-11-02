#include "pattern_index.h"
#include "representation.h"
#include "common/cli.h"
#include "common/index_library.h"
#include "common/log.h"
#include <chrono>

using namespace mcviz;
using Clock = std::chrono::steady_clock;

int main(int argc, char** argv) {
    const std::string mode = argv[1];
    const auto engine = mode == "exact" ? Engine::exact : Engine::upper;
    CliCapabilities caps;
    caps.commands = {"query"};
    const auto opt = parse_cli(argc - 1, argv + 1, caps);
    auto g = Graph::load(opt.layers, opt.graph_file, opt.vertices_file);
    auto path = opt.index_file;
    if (path.empty()) {
        const auto spec = build_spec(opt);
        const auto paths = library_paths(opt, g, spec, mode == "exact" ? "MCviz" : "MCviz-plus");
        path = find_library_index(paths, g, spec, engine_name(engine));
    }
    auto idx = PatternIndex::load(g, path, engine);
    const auto spec = idx.identity.spec;

    const auto tree = prepare_representation(g, idx.forest);
    g = Graph{};
    idx.weights = {};
    idx.weights.shrink_to_fit();

    idx.core_domains = {};

    idx.forest = Forest{};

    auto start = Clock::now();
    const auto linear = tree_to_linear(tree);
    const double forward = std::chrono::duration<double>(Clock::now() - start).count();
    start = Clock::now();
    const auto restored = linear_to_tree(linear);
    const double reverse = std::chrono::duration<double>(Clock::now() - start).count();
    log_conversion(opt, spec, mode == "exact" ? "mcviz" : "mcviz-plus",
                   forward, reverse, representation_bytes(tree), representation_bytes(linear));
}
