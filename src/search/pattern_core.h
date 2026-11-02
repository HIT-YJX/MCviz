#pragma once
#include "pattern_bounds.h"

namespace mcviz {
enum class CoreBackend { reference, queue };
struct CoreMetrics {
    std::uint64_t pair_probes = 0, adjacency_scans = 0, degree_updates = 0, queue_pops = 0;
    std::uint64_t slots = 0, arcs = 0, map_initializations = 0;
};

class CoreWorkspace {
public:
    void bind(const Graph& g, const Sets& domains, CoreMetrics* metrics = nullptr);
    std::optional<std::size_t> find(Layer i, VertexPos u) const;
private:
    std::vector<std::size_t> offsets_, ids_;
    std::vector<std::uint64_t> stamps_;
    std::uint64_t epoch_ = 0;
};
struct CoreResult {
    Sets domains;
    bool possible = true;
    std::uint64_t removed = 0, rounds = 0;
};
CoreResult reduce(const Graph& g, const LayerPattern& p, const Thresholds& gamma,
                  Sets domains, const Required& required = {},
                  CoreBackend backend = CoreBackend::reference, CoreWorkspace* workspace = nullptr,
                  CoreMetrics* metrics = nullptr);
CoreResult reduce_queue(const Graph& g, const LayerPattern& p, const Thresholds& gamma,
                        Sets domains, const Required& required, CoreWorkspace* workspace,
                        CoreMetrics* metrics);
}
