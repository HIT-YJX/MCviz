#pragma once
#include "contract.h"
#include <stdexcept>

namespace mcviz {

struct CliCapabilities {
    std::vector<std::string> commands;
    bool solver = false, edge_bound = false, method = false, min_layers = false;
    bool cohesion_build = false;
    bool support_build = false;
};


struct CliOptions {
    std::string command, graph_file, vertices_file, index_file;
    std::string index_dir, log_dir;
    unsigned layers = 0;
    std::optional<std::string> pattern_selector;
    std::string path_order_file;
    std::optional<Layer> primary;
    std::optional<Objective> objective;
    std::optional<Thresholds> beta;
    QueryRequest query;
    std::string solver = "auto", edge_bound = "auto", method;
    std::optional<unsigned> min_layers;
    unsigned build_threads = 1;
    bool build_threads_explicit = false;
    std::uint64_t query_heuristic_ms = 1000;
    std::uint64_t pivot_memory_bytes = 64ULL << 20;
    std::uint64_t cohesion_work_budget = 20000000;
    std::uint64_t cohesion_cache_bytes = 64ULL << 20;
    std::string cohesion_strategy = "auto", build_filter = "auto";
    bool help = false;
};

CliOptions parse_cli(int argc, const char* const* argv, const CliCapabilities& capabilities);
BuildSpec build_spec(const CliOptions& options);
std::string cli_help();

}
