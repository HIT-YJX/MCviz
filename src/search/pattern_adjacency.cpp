#include "pattern_adjacency.h"
#include "common/text.h"
#include <algorithm>
#include <chrono>
#include <limits>
#include <new>
#include <stdexcept>

namespace mcviz {
namespace {
using Clock = std::chrono::steady_clock;
constexpr auto missing = std::numeric_limits<std::size_t>::max();


// Check complete logical sizes before allocation; overflow always refuses a
// cache, even when the caller supplied the maximum integer budget. O(1) work.
struct Size {
    std::uint64_t bytes = 0;
    bool overflow = false;
    void add(std::uint64_t count, std::uint64_t width) {
        if (count && width > (UINT64_MAX - bytes) / count) {
            overflow = true;
            return;
        }
        bytes += count * width;
    }
    bool fits(std::uint64_t limit) const { return !overflow && bytes <= limit; }
};

// Each target layer has this many words; Graph positions fit unsigned.
// Time and auxiliary space O(1).
std::size_t words(std::size_t count) { return count / 64 + (count % 64 != 0); }
} // namespace

// Build one row per base vertex for every directed P edge. Let N be the full
// graph size, n_i=|base_i|, b_i=ceil(n_i/64), h=sum_i n_i*deg_P(i), and a the
// scanned original adjacency entries. Time O(N+K+|P|+sum_i n_i log(n_i+1)
// +h log(K+1)+a+S), space
// O(N+K+|P|+sum_i n_i+S), S=sum_(i,j directed P) n_i*b_j. No subset states,
// bipartition, connectivity or star assumptions enter the construction.
PatternAdjacency::PatternAdjacency(const Graph& graph, const LayerPattern& pattern,
                                   const Sets& base, bool enabled, std::uint64_t cache_bytes)
    : graph_(graph), pattern_(pattern), base_(base) {
    const auto start = Clock::now();
    // Borrow an already ordered base; only noncanonical caller input needs a
    // normalized owned copy. Validation costs O(n) once, never once per edge.
    for (const auto& row : base) if (!std::is_sorted(row.begin(), row.end())) {
        ordered_ = base;
        for (auto& a : ordered_) std::sort(a.begin(), a.end());
        owns_sparse_base_ = true;
        break;
    }
    for (Layer i = 0; i < graph.K; ++i) {
        const auto& a = owns_sparse_base_ ? ordered_[i] : base[i];
        for (auto u : a) graph.get_id(i, u);
    }
    auto prepare = [&] {

        if (!enabled) return;
        Size size;
        size.add(1, sizeof(PatternAdjacency));
        size.add(graph.K, sizeof(std::vector<VertexPos>));
        size.add(graph.K, 2 * sizeof(std::vector<std::size_t>));
        std::size_t count = 0;
        for (Layer i = 0; i < graph.K; ++i) {
            size.add(base[i].size(), sizeof(VertexPos));
            size.add(graph.vertex_ids[i].size(), sizeof(std::size_t));
            size.add(pattern.neighbors[i].size(), sizeof(std::size_t) + sizeof(Block));
            if (pattern.neighbors[i].size() > missing - count) return;
            count += pattern.neighbors[i].size();
            for (auto j : pattern.neighbors[i]) {
                size.add(base[i].size(), words(base[j].size()) * sizeof(std::uint64_t));
            }
        }
        if (!size.fits(cache_bytes)) return;

        if (!owns_sparse_base_) ordered_ = base;
        position_.resize(graph.K);
        outgoing_.resize(graph.K);
        blocks_.reserve(count);
        for (Layer i = 0; i < graph.K; ++i) {
            auto& layer = ordered_[i];
            std::sort(layer.begin(), layer.end());
            position_[i].assign(graph.vertex_ids[i].size(), missing);
            for (std::size_t r = 0; r < layer.size(); ++r) {
                graph.get_id(i, layer[r]);
                position_[i][layer[r]] = r;
            }
            outgoing_[i].reserve(pattern.neighbors[i].size());
            for (auto j : pattern.neighbors[i]) {
                Block block;
                block.from = i;
                block.to = j;
                block.words = words(base[j].size());
                if (block.words && layer.size() > missing / block.words) return;
                block.bits.assign(layer.size() * block.words, 0);
                outgoing_[i].push_back(blocks_.size());
                blocks_.push_back(std::move(block));
            }
        }

        // Count actual owned capacities, including all nested vector headers.
        // This catches capacity rounding before committing the ready flag.
        Size actual;
        actual.add(1, sizeof(PatternAdjacency));
        actual.add(ordered_.capacity(), sizeof(std::vector<VertexPos>));
        actual.add(position_.capacity(), sizeof(std::vector<std::size_t>));
        actual.add(outgoing_.capacity(), sizeof(std::vector<std::size_t>));
        actual.add(blocks_.capacity(), sizeof(Block));
        for (Layer i = 0; i < graph.K; ++i) {
            actual.add(ordered_[i].capacity(), sizeof(VertexPos));
            actual.add(position_[i].capacity(), sizeof(std::size_t));
            actual.add(outgoing_[i].capacity(), sizeof(std::size_t));
        }
        for (const auto& block : blocks_) actual.add(block.bits.capacity(), sizeof(std::uint64_t));
        if (!actual.fits(cache_bytes)) return;

        // Original adjacency determines the bits; missing non-P edges are never
        // queried. All maps are installed before any directed block is filled.
        for (auto& block : blocks_) {
            if (!block.words) continue;
            const auto& from = ordered_[block.from];
            const auto& to = position_[block.to];
            for (std::size_t r = 0; r < from.size(); ++r)
                for (auto v : graph.neighbors(block.from, from[r], block.to)) {
                    const auto at = to[v];
                    if (at != missing) block.bits[r * block.words + at / 64] |= 1ULL << (at % 64);
                }
        }
        logical_bytes_ = actual.bytes;
        ready_ = true;
    
    };
    prepare();
    if (!ready_) {
        if (!owns_sparse_base_) Sets().swap(ordered_);
        std::vector<std::vector<std::size_t>>().swap(position_);
        std::vector<std::vector<std::size_t>>().swap(outgoing_);
        std::vector<Block>().swap(blocks_);
        logical_bytes_ = 0;
        // A noncanonical input's normalization copy is compulsory sparse
        // storage, not an admitted bit cache; expose it even at budget zero.
        if (owns_sparse_base_) {
            logical_bytes_ = ordered_.capacity() * sizeof(std::vector<VertexPos>);
            for (const auto& a : ordered_) logical_bytes_ += a.capacity() * sizeof(VertexPos);
        }
    }
    prepare_ms_ = std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

// Object identity prevents a cache from serving another Graph/P with equal sizes.
// Compare the complete normalized base, including isolated layers. O(K+n) time
// with a ready cache or sorted borrowed base; otherwise O(sum_i n_i log(n_i+1))
// time and O(K+n) temporary storage. No cache/input state is changed.
bool PatternAdjacency::matches(const Graph& graph, const LayerPattern& pattern,
                               const Sets& base) const {
    if (&graph != &graph_ || &pattern != &pattern_) return false;
    if (ready_) return ordered_ == base;
    return (owns_sparse_base_ ? ordered_ : base_) == base;
}

// Intersect the shortest of base[j] and all mandatory neighbor lists. Sorted
// lists provide exact membership without a graph-sized per-edge scan or map.
// With r_j constrained lists and shortest length s_j, time is O(K+|Q|K+
// sum_j s_j*(log(n_j+1)+sum_q log(deg_j(q)+1))), plus untouched-layer output;
// scratch is O(|Q|), and explicit output includes every isolated layer.
DomainResult PatternAdjacency::condition_sparse(const Required& required) const {
    for (auto q : required) graph_.get_id(q.layer, q.pos);
    DomainResult out{Sets(graph_.K), true};
    const auto& base = owns_sparse_base_ ? ordered_ : base_;
    for (Layer j = 0; j < graph_.K; ++j) {
        std::vector<const std::vector<VertexPos>*> rows{&base[j]};
        for (auto q : required) if (pattern_.adjacent(q.layer, j))
            rows.push_back(&graph_.neighbors(q.layer, q.pos, j));
        if (rows.size() == 1) { out.domains[j] = base[j]; continue; }
        const auto seed = std::min_element(rows.begin(), rows.end(),
            [](auto a, auto b) { return a->size() < b->size(); });
        std::swap(rows.front(), *seed);
        auto& d = out.domains[j];
        d.reserve(rows.front()->size());
        for (auto u : *rows.front()) {
            bool keep = true;
            for (std::size_t r = 1; r < rows.size() && keep; ++r)
                keep = std::binary_search(rows[r]->begin(), rows[r]->end(), u);
            if (keep) d.push_back(u);
        }
    }
    out.possible = contains_required(out.domains, required);
    return out;
}

// Each mandatory point contributes one AND in each of its P-neighbor layers.
// A shared target of multiple mandatory layers is intersected multiple times,
// so triangles/cycles obey the same conjunction as general condition_domains.
// With cache: O(K+sum_i n_i+sum_(q,j in N_P(q.layer)) b_j+|Q| log(n+1)) time,
// output O(K+sum_i n_i), scratch O(K+sum_touched_j b_j). Otherwise sparse
// intersection also handles graph-valid mandatory points outside the base.
DomainResult PatternAdjacency::condition(const Required& required) const {
    if (!ready_) return condition_sparse(required);
    bool inside = true;
    for (auto q : required) {
        graph_.get_id(q.layer, q.pos);
        if (position_[q.layer][q.pos] == missing) inside = false;
    }
    if (!inside) return condition_sparse(required);

    DomainResult out{Sets(graph_.K), true};
    std::vector<std::vector<std::uint64_t>> masks(graph_.K);
    std::vector<unsigned char> touched(graph_.K, 0);
    for (auto q : required) {
        const auto row = position_[q.layer][q.pos];
        for (auto id : outgoing_[q.layer]) {
            const auto& block = blocks_[id];
            auto& mask = masks[block.to];
            if (!touched[block.to]) {
                mask.assign(block.words, UINT64_MAX);
                const auto tail = ordered_[block.to].size() % 64;
                if (tail) mask.back() = (1ULL << tail) - 1;
                touched[block.to] = 1;
            }
            for (std::size_t k = 0; k < block.words; ++k)
                mask[k] &= block.bits[row * block.words + k];
        }
    }
    for (Layer j = 0; j < graph_.K; ++j) {
        if (!touched[j]) { out.domains[j] = ordered_[j]; continue; }
        auto& domain = out.domains[j];
        domain.clear();
        for (std::size_t k = 0; k < masks[j].size(); ++k) {
            auto bits = masks[j][k];
            while (bits) {
                const auto r = 64 * k + static_cast<unsigned>(__builtin_ctzll(bits));
                domain.push_back(ordered_[j][r]);
                bits &= bits - 1;
            }
        }
    }
    out.possible = contains_required(out.domains, required);
    return out;
}
} // namespace mcviz
