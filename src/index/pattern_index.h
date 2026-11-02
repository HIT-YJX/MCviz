#pragma once
#include "forest.h"
#include "search/pattern_solver.h"
#include "common/output.h"
#include <stdexcept>

namespace mcviz {
enum class Engine { exact, upper };
struct IndexOptions {
    std::string solver = "auto", edge_bound = "auto";
    SolverOptions search = [] { SolverOptions opt; opt.pivot_memory_bytes = 512ULL << 20; return opt; }();
    std::uint64_t cohesion_work_budget = 20000000, cohesion_cache_bytes = 64ULL << 20;
    CoreBackend build_core_backend = CoreBackend::queue;
    std::string build_filter = "auto";
    std::uint64_t build_filter_work_budget = 100000000, build_filter_cache_bytes = 64ULL << 20;
    bool build_capacity_degrees = true;
    bool optimize_star_build = true;
    unsigned build_threads = 1;
    std::uint64_t support_cache_bytes = 64ULL << 20;
    std::uint64_t query_heuristic_ms = 1000;
};
class PatternIndex {
public:
    Engine engine = Engine::exact;
    BuildIdentity identity;
    IndexOptions options;
    std::string execution_path, bound;
    std::string source_graph_sha256, source_vertices_sha256;
    std::vector<WeightedEdge> weights;
    Forest forest;
    Sets core_domains;
    SolverStats build_stats;
    std::uint64_t build_solver_calls = 0;
    std::string build_strategy;
    std::map<std::string, double> build_timing;
    std::map<std::string, std::uint64_t> build_counters;

    static PatternIndex construct(const Graph& g, const BuildSpec& spec, Engine engine,
                                  const IndexOptions& options = {});
    SearchResult query(const Graph& g, const QueryRequest& q, RunInfo& run) const;
    std::string query_execution_path() const;
    void save(const Graph& g, const std::string& path) const;
    static PatternIndex load(const Graph& g, const std::string& path, Engine engine);
};
const char* engine_name(Engine engine);
}
