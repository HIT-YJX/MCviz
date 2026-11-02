#include "star_support_build.h"
#include "common/text.h"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <thread>

namespace mcviz {
namespace {
using Clock = std::chrono::steady_clock;
using NeighborList = std::vector<VertexPos>;

std::uint64_t add_count(std::uint64_t a, std::uint64_t b) {
    return b > std::numeric_limits<std::uint64_t>::max() - a
        ? std::numeric_limits<std::uint64_t>::max() : a + b;
}
std::uint64_t multiply_count(std::uint64_t a, std::uint64_t b) {
    return a && b > std::numeric_limits<std::uint64_t>::max() / a
        ? std::numeric_limits<std::uint64_t>::max() : a * b;
}
Score add_score(Score a, Score b) {
    if (b > std::numeric_limits<Score>::max() - a) throw std::runtime_error("Score addition overflow");
    return a + b;
}
Score multiply_score(Score a, Score b) {
    if (a && b > std::numeric_limits<Score>::max() / a) throw std::runtime_error("Score multiplication overflow");
    return a * b;
}
double milliseconds(Clock::time_point start) {
    return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}
Layer center_of(const LayerPattern& pattern) {
    auto centers = pattern.star_centers();

    return std::find(centers.begin(), centers.end(), pattern.anchor_layer) == centers.end()
        ? centers.front() : pattern.anchor_layer;
}

struct Incident { VertexPos primary; std::size_t edge; };
struct Row {
    Layer layer = 0;
    VertexPos vertex = 0;
    const NeighborList* primary = nullptr;
    std::vector<Incident> incident;
    std::uint64_t cost = 0;
};

struct Topology {
    Layer center;
    unsigned layers;
    std::vector<std::size_t> offsets, order, starts;
    std::vector<Row> rows;
    std::vector<std::vector<std::size_t>> primary_neighbors;

    Topology(const Graph& graph, Layer primary_layer, std::size_t secondary_count)
        : center(primary_layer), layers(graph.K), offsets(layers), rows(secondary_count),
          primary_neighbors(graph.vertex_ids[center].size()) {
        std::size_t next = 0;
        for (Layer j = 0; j < layers; ++j) {
            offsets[j] = next;
            if (j == center) continue;
            for (VertexPos v = 0; v < graph.vertex_ids[j].size(); ++v) {
                auto& row = rows[next++];
                row.layer = j;
                row.vertex = v;
                row.primary = &graph.neighbors(j, v, center);
            }
        }

        for (std::size_t id = 0; id < graph.edges.size(); ++id) {
            auto edge = graph.edges[id];
            if (edge.u.layer != center) std::swap(edge.u, edge.v);

            const auto flat = offsets[edge.v.layer] + edge.v.pos;
            primary_neighbors[edge.u.pos].push_back(flat);
            rows[flat].incident.push_back({edge.u.pos, id});
        }
        require(primary_neighbors.size() <= std::numeric_limits<std::size_t>::max() / (std::size_t(layers) + 1),
                "Support primary range size overflow");
        starts.resize(primary_neighbors.size() * (std::size_t(layers) + 1));
        for (std::size_t u = 0; u < primary_neighbors.size(); ++u) {
            auto& neighbors = primary_neighbors[u];
            std::sort(neighbors.begin(), neighbors.end());
            std::size_t at = 0;
            for (Layer j = 0; j < layers; ++j) {
                starts[u * (std::size_t(layers) + 1) + j] = at;
                while (at < neighbors.size() && rows[neighbors[at]].layer == j) ++at;
            }
            starts[u * (std::size_t(layers) + 1) + layers] = at;
        }
        for (auto& row : rows) {
            for (auto u : *row.primary) row.cost = add_count(row.cost, primary_neighbors[u].size());
        }
        order.resize(rows.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](auto a, auto b) {
            return rows[a].cost != rows[b].cost ? rows[a].cost > rows[b].cost : a < b;
        });
    }

    std::uint64_t bytes() const {
        auto result = multiply_count(offsets.capacity() + order.capacity() + starts.capacity(), sizeof(std::size_t));
        result = add_count(result, multiply_count(rows.capacity(), sizeof(Row)));
        result = add_count(result, multiply_count(primary_neighbors.capacity(), sizeof(std::vector<std::size_t>)));
        for (const auto& row : rows) result = add_count(result, multiply_count(row.incident.capacity(), sizeof(Incident)));
        for (const auto& neighbors : primary_neighbors) result = add_count(result, multiply_count(neighbors.capacity(), sizeof(std::size_t)));
        return result;
    }
};

struct Worker {
    std::unique_ptr<VertexPos[]> common;
    std::unique_ptr<std::size_t[]> touched;
    std::unique_ptr<Score[]> histogram, sums, minima, capacity;
    std::unique_ptr<unsigned char[]> feasible;
    std::size_t touched_size = 0;
    std::uint64_t payload_bytes;
    std::uint64_t rows = 0, edges = 0, wedges = 0, merges = 0, lookups = 0;
    double rows_ms = 0, weights_ms = 0, total_ms = 0;

    Worker(std::size_t secondary, std::size_t maximum_degree, unsigned layers)
        : common(new VertexPos[secondary]()),
          touched(new std::size_t[secondary]),
          histogram(new Score[maximum_degree + 1]), sums(new Score[maximum_degree + 1]),
          minima(new Score[maximum_degree + 1]), capacity(new Score[layers]),
          feasible(new unsigned char[maximum_degree + 1]),
          payload_bytes(estimate(secondary, maximum_degree, layers)) {}
    static std::uint64_t estimate(std::size_t secondary, std::size_t degree, unsigned layers) {
        auto bytes = multiply_count(degree + 1, 3 * sizeof(Score) + sizeof(unsigned char));
        bytes = add_count(bytes, multiply_count(layers, sizeof(Score)));
        bytes = add_count(bytes, multiply_count(secondary, sizeof(VertexPos) + sizeof(std::size_t)));
        return bytes;
    }
    std::uint64_t bytes() const { return payload_bytes; }
};

// Count shared primary neighbors once for every incident edge of this secondary row.
// Time O(t + sum_{u in N(row)} deg(u)), where t is the previous touched count; O(n) cache.
void prepare_row(const Topology& topology, const Row& row, Worker& worker) {
    for (std::size_t i = 0; i < worker.touched_size; ++i) worker.common[worker.touched[i]] = 0;
    worker.touched_size = 0;
    for (auto u : *row.primary) for (auto x : topology.primary_neighbors[u]) {
        if (!worker.common[x]) worker.touched[worker.touched_size++] = x;
        ++worker.common[x];
    }
    ++worker.rows;
    worker.wedges = add_count(worker.wedges, row.cost);
}

Score evaluate(const Topology& topology, const Row& row, VertexPos u, const LayerPattern& pattern,
               Objective objective, const Thresholds& beta, Worker& worker) {
    const auto maximum = row.primary->size();
    const auto base = std::size_t(u) * (std::size_t(topology.layers) + 1);
    const auto& neighbors = topology.primary_neighbors[u];
    worker.capacity[topology.center] = maximum;
    for (Layer j = 0; j < topology.layers; ++j) if (j != topology.center)
        worker.capacity[j] = topology.starts[base + j + 1] - topology.starts[base + j];
    for (Layer j = 0; j < topology.layers; ++j) if (worker.capacity[j] < beta[j]) return 0;
    Score coarse = 0;
    if (objective == Objective::edge) {
        for (auto edge : pattern.edges)
            coarse = add_score(coarse, multiply_score(worker.capacity[edge.first], worker.capacity[edge.second]));
    } else if (objective == Objective::vertex) {
        for (Layer j = 0; j < topology.layers; ++j) coarse = add_score(coarse, worker.capacity[j]);
    } else coarse = *std::min_element(worker.capacity.get(), worker.capacity.get() + topology.layers);
    std::fill_n(worker.sums.get(), maximum + 1, 0);
    std::fill_n(worker.minima.get(), maximum + 1, std::numeric_limits<Score>::max());
    std::fill_n(worker.feasible.get(), maximum + 1, 1);
    for (Layer j = 0; j < topology.layers; ++j) if (j != topology.center) {
        std::fill_n(worker.histogram.get(), maximum + 1, 0);
        for (std::size_t i = topology.starts[base + j]; i < topology.starts[base + j + 1]; ++i) {
            const auto x = neighbors[i];
            const std::size_t count = worker.common[x] - 1;
            ++worker.histogram[count];
        }
        Score count = 0;
        for (std::size_t a = maximum; a > 0; --a) {
            count += worker.histogram[a - 1];
            worker.sums[a] = add_score(worker.sums[a], count);
            worker.minima[a] = std::min(worker.minima[a], count);
            if (count < beta[j]) worker.feasible[a] = 0;
        }
    }
    Score bound = 0;
    for (std::size_t a = static_cast<std::size_t>(beta[topology.center]); a <= maximum; ++a)
        if (worker.feasible[a]) {
            const auto value = objective == Objective::edge ? multiply_score(a, worker.sums[a]) :
                objective == Objective::vertex ? add_score(a, worker.sums[a]) : std::min<Score>(a, worker.minima[a]);
            bound = std::max(bound, value);
        }
    return std::min(coarse, bound);
}

}

// Reuse typed support rows to bound each edge, with scalar fallback when cache admission fails.
// Setup stores O(n+m+K) topology; work sums row intersections and support histograms.
// Each worker owns its bounded cache and output workspace; budget refusal uses scalar bounds.
SupportBuildResult build_star_support(const Graph& active, const LayerPattern& pattern,
                                     Objective objective, const Thresholds& beta,
                                     unsigned threads, std::uint64_t cache_bytes) {
    SupportBuildResult result;
    result.values.resize(active.edges.size());
    for (auto key : {"support_requested_threads", "support_actual_threads", "support_workers_allocated",
                     "support_rows_prepared", "support_edges_evaluated", "support_cache_peak_bytes",
                     "support_cache_budget_bytes", "support_scratch_bytes_per_worker", "support_shared_bytes",
                     "support_fallback_rows", "support_thread_fallbacks", "support_memory_fallbacks",
                     "support_wedge_visits", "support_merge_pairs", "support_pair_lookups", "support_scalar_edges"})
        result.counters[key] = 0;
    for (auto key : {"support_prepare_ms", "support_rows_ms", "support_weights_ms"}) result.timing[key] = 0;
    result.counters["support_requested_threads"] = threads;
    result.counters["support_cache_budget_bytes"] = cache_bytes;
    auto start = Clock::now();
    PreparedStarBounds scalar(active, pattern, objective, beta);
    const auto center = center_of(pattern);
    std::size_t secondary_count = 0, maximum_degree = 0;
    for (Layer j = 0; j < active.K; ++j) if (j != center) {
        require(active.vertex_ids[j].size() <= std::numeric_limits<std::size_t>::max() - secondary_count,
                "Support secondary size overflow");
        secondary_count += active.vertex_ids[j].size();
        for (VertexPos v = 0; v < active.vertex_ids[j].size(); ++v)
            maximum_degree = std::max(maximum_degree, active.neighbors(j, v, center).size());
    }
    require(maximum_degree < std::numeric_limits<std::size_t>::max(), "Support histogram size overflow");
    const auto estimate = Worker::estimate(secondary_count, maximum_degree, active.K);
    result.counters["support_scratch_bytes_per_worker"] = estimate;
    auto scalar_fallback = [&]() {
        result.strategy = "star-support-scalar";
        result.counters["support_actual_threads"] = 1;
        result.counters["support_fallback_rows"] = secondary_count;
        result.counters["support_scalar_edges"] = active.edges.size();
        result.counters["support_edges_evaluated"] = active.edges.size();
        result.timing["support_prepare_ms"] = milliseconds(start);
        const auto weights_start = Clock::now();
        for (std::size_t id = 0; id < active.edges.size(); ++id)
            result.values[id] = scalar.bound(active.edges[id], true);
        result.timing["support_weights_ms"] = milliseconds(weights_start);
        result.timing["support_parallel_wall_ms"] = result.timing["support_weights_ms"];
        result.timing["worker_0_total_ms"] = result.timing["support_weights_ms"];
        result.timing["worker_0_prepare_ms"] = 0;
        result.timing["worker_0_solver_ms"] = result.timing["support_weights_ms"];
        result.counters["worker_0_tasks"] = active.edges.size();
        result.counters["worker_0_logical_bytes"] = 0;
    };
    if (estimate > cache_bytes || !secondary_count) {
        scalar_fallback();
        return result;
    }
    std::unique_ptr<Topology> topology;
    std::vector<std::unique_ptr<Worker>> workers;
    auto count = static_cast<std::size_t>(std::min<std::uint64_t>({threads, secondary_count, cache_bytes / estimate}));
    std::uint64_t allocated = 0, reserved_peak = 0;
    {
        topology = std::make_unique<Topology>(active, center, secondary_count);
        workers.reserve(count);
        for (std::size_t i = 0; i < count; ++i) {
            reserved_peak = std::max(reserved_peak, add_count(allocated, estimate));
            auto worker = std::make_unique<Worker>(secondary_count, maximum_degree, active.K);
            const auto bytes = worker->bytes();
            if (bytes > cache_bytes - allocated) break;
            allocated += bytes;
            result.counters["support_scratch_bytes_per_worker"] = std::max(estimate, bytes);
            workers.push_back(std::move(worker));
        }
    }
    result.counters["support_cache_peak_bytes"] = reserved_peak;
    if (workers.empty()) {
        topology.reset();
        scalar_fallback();
        return result;
    }
    result.strategy = "star-support-row-cache";
    result.counters["support_workers_allocated"] = workers.size();
    result.counters["support_shared_bytes"] = topology->bytes();
    result.timing["support_prepare_ms"] = milliseconds(start);
    std::atomic<std::size_t> next{0};

    auto run = [&](std::size_t worker_id) {
        auto& worker = *workers[worker_id];
        const auto worker_start = Clock::now();
        for (;;) {
            const auto task = next.fetch_add(1, std::memory_order_relaxed);
            if (task >= topology->order.size()) break;
            const auto& row = topology->rows[topology->order[task]];
            const auto row_start = Clock::now();
            prepare_row(*topology, row, worker);
            worker.rows_ms += milliseconds(row_start);
            worker.lookups = add_count(worker.lookups, row.cost);
            const auto weights_start = Clock::now();
            for (const auto edge : row.incident) {
                result.values[edge.edge] = evaluate(*topology, row, edge.primary, pattern, objective, beta, worker);
                ++worker.edges;
            }
            worker.weights_ms += milliseconds(weights_start);
        }
        worker.total_ms = milliseconds(worker_start);
    };
    std::vector<std::thread> pool;
    std::uint64_t thread_fallbacks = 0;
    const auto parallel_start = Clock::now();
    { pool.reserve(workers.size() - 1); }
    for (std::size_t i = 0; i + 1 < workers.size(); ++i) {
        { pool.emplace_back(run, i); }
    }
    run(pool.size());
    for (auto& thread : pool) if (thread.joinable()) thread.join();
    result.timing["support_parallel_wall_ms"] = milliseconds(parallel_start);
    result.counters["support_actual_threads"] = pool.size() + 1;
    result.counters["support_thread_fallbacks"] = thread_fallbacks;
    for (const auto& worker : workers) {
        result.counters["support_rows_prepared"] = add_count(result.counters["support_rows_prepared"], worker->rows);
        result.counters["support_edges_evaluated"] = add_count(result.counters["support_edges_evaluated"], worker->edges);
        result.counters["support_wedge_visits"] = add_count(result.counters["support_wedge_visits"], worker->wedges);
        result.counters["support_merge_pairs"] = add_count(result.counters["support_merge_pairs"], worker->merges);
        result.counters["support_pair_lookups"] = add_count(result.counters["support_pair_lookups"], worker->lookups);
        result.timing["support_rows_ms"] += worker->rows_ms;
        result.timing["support_weights_ms"] += worker->weights_ms;
    }
    for (std::size_t i = 0; i <= pool.size(); ++i) {
        const auto key = "worker_" + std::to_string(i);
        const auto& worker = *workers[i];
        result.timing[key + "_total_ms"] = worker.total_ms;
        result.timing[key + "_prepare_ms"] = worker.rows_ms;
        result.timing[key + "_solver_ms"] = worker.weights_ms;
        result.counters[key + "_tasks"] = worker.edges;
        result.counters[key + "_logical_bytes"] = worker.bytes();
    }
    return result;
}
}
