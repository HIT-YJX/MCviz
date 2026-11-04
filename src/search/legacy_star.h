#pragma once

#include "common/graph.h"
#include "common/objective.h"
#include <algorithm>
#include <climits>
#include <deque>
#include <functional>
#include <limits>
#include <numeric>
#include <stdexcept>
namespace mcviz::legacy_star {
using namespace std;
using uint = unsigned;
using ull = unsigned long long;
using KPartiteGraph = mcviz::Graph;
enum ProblemType { MESC = 1, MVSC = 2, MBSC = 3 };
struct StarClique {
    vector<vector<uint>> sets;
    ull size = 0;
    explicit StarClique(uint K = 0) : sets(K) {}
};

inline ull compute_star_size(const vector<uint>& a, ProblemType obj, uint p = 0) {
    if (obj != MESC && obj != MVSC && obj != MBSC) throw invalid_argument("Invalid star objective");
    if (a.empty()) return 0;
    if (p >= a.size()) throw invalid_argument("Invalid star primary");
    ull sum = 0;
    for (uint i = 0; i < a.size(); ++i)
        if (i != p) sum = mcviz::checked_add(sum, a[i]);
    if (obj == MESC) return mcviz::checked_mul(a[p], sum);
    if (obj == MVSC) return mcviz::checked_add(a[p], sum);
    return *min_element(a.begin(), a.end());
}
// Enumerate star witnesses with two-sided pivots and shared budget accounting.
template <class Charge> class ExactPivotEnumeration {
    using Bits = vector<ull>;
    const KPartiteGraph& g;
    uint p;
    ProblemType objective;
    const vector<uint>& beta;
    const vector<uint>& primary;
    const vector<pair<uint, uint>>& secondary;
    const vector<Bits>& rows;
    const vector<Bits>& columns;
    StarClique& best;
    const function<void(const StarClique&)>& publish;
    const function<bool(ull, const vector<pair<uint, uint>>&, const vector<uint>&)>& improve;
    ull& nodes;
    Charge charge;
    bool& refused;
    size_t nu, nv, wu, wv;
    vector<uint> selected_u, selected_v, selected_sizes;
    Bits selected_u_bits;
    vector<vector<pair<size_t, ull>>> layer_words;
    struct Frame {
        Bits pu, pv, xu, xv, allowed;
        vector<uint> right, sizes, hist, branch_u, branch_v, support;
        vector<uint> refine_layers;
        vector<pair<size_t, ull>> refine_words;
        vector<size_t> refine_offsets;
        vector<unsigned char> refine_dirty;
        bool refine_support_ready = false;
    };
    deque<Frame> frames;

    bool prepared = false;
    size_t maximum_primary = 0;

    static void put(Bits& bits, size_t v) { bits[v / 64] |= 1ULL << (v % 64); }

    static void erase(Bits& bits, size_t v) { bits[v / 64] &= ~(1ULL << (v % 64)); }



    bool empty(const Bits& bits) {
        for (ull word : bits) {
            if (!charge(1)) {
                return false;
            }
            if (word) return false;
        }
        return true;
    }

    uint count(const Bits& bits) {
        uint result = 0;
        for (ull word : bits) {
            if (!charge(1)) {
                return 0;
            }
            result += __builtin_popcountll(word);
        }
        return result;
    }


    template <class F> void each(const Bits& bits, F fn) {
        for (size_t w = 0; w < bits.size(); ++w) {
            if (!charge(1)) {
                return;
            }
            ull word = bits[w];
            while (word) {
                uint bit = __builtin_ctzll(word);
                word &= word - 1;
                if (!charge(1)) {
                    return;
                }
                fn(static_cast<uint>(64 * w + bit));
                if (refused) return;
            }
        }
    }



    void copy(Bits& out, const Bits& in) {
        out.resize(in.size());
        for (size_t w = 0; w < in.size(); ++w) {
            if (!charge(1)) {
                return;
            }
            out[w] = in[w];
        }
    }



    void intersect(Bits& out, const Bits& in, const Bits& mask) {
        out.resize(in.size());
        for (size_t w = 0; w < in.size(); ++w) {
            if (!charge(1)) {
                return;
            }
            out[w] = in[w] & mask[w];
        }
    }



    uint intersection_count(const Bits& a, const Bits& b) {
        uint result = 0;
        for (size_t w = 0; w < a.size(); ++w) {
            if (!charge(1)) {
                return 0;
            }
            result += __builtin_popcountll(a[w] & b[w]);
        }
        return result;
    }


    void fill_all(Bits& bits, size_t vertices) {
        bits.resize((vertices + 63) / 64);
        for (auto& word : bits) {
            if (!charge(1)) {
                return;
            }
            word = ~0ULL;
        }
        if (vertices % 64) bits.back() = (1ULL << (vertices % 64)) - 1;
    }

    void initialize(Frame& f) {
        fill_all(f.pu, nu);
        fill_all(f.pv, nv);
        f.xu.assign(wu, 0);
        f.xv.assign(wv, 0);
        if (!charge(wu + wv)) {
            return;
        }
    }



    // Peel candidates that cannot meet a neighboring layer quota.
    bool refine_compact(Frame& f, uint& np, uint& ns) {
        f.refine_layers.clear();
        f.refine_dirty.assign(g.K, 0);
        if (!charge(g.K)) {
            return false;
        }
        for (uint j = 0; j < g.K; ++j) {
            if (!charge(1)) {
                return false;
            }
            if (j != p && selected_sizes[j] < beta[j]) {
                f.refine_layers.push_back(j);
                f.refine_dirty[j] = 1;
            }
        }
        const bool check_columns = selected_sizes[p] < beta[p];
        if (f.refine_layers.empty() && !check_columns) return true;
        bool changed = false;
        auto discard_column = [&](uint v) {
            erase(f.pv, v);
            f.refine_dirty[secondary[v].first] = 1;
            changed = true;
        };

        if (check_columns) {
            if (f.refine_support_ready) {
                for (size_t i = selected_v.size(); i < f.right.size(); ++i) {
                    if (!charge(1)) {
                        return false;
                    }
                    if (f.support[i] < beta[p]) discard_column(f.right[i]);
                }
            } else
                each(f.pv, [&](uint v) {
                    uint degree = selected_sizes[p];
                    for (size_t w = 0; w < wu && degree < beta[p]; ++w) {
                        if (!charge(1)) {
                            return;
                        }
                        degree += __builtin_popcountll(columns[v][w] & f.pu[w]);
                    }
                    if (degree < beta[p]) discard_column(v);
                });
        }
        if (refused) return false;
        f.refine_support_ready = false;
        while (!f.refine_layers.empty()) {

            f.refine_words.clear();
            f.refine_offsets.assign(g.K + 1, 0);
            if (!charge(g.K + 1)) {
                return false;
            }
            for (uint j : f.refine_layers) {
                f.refine_offsets[j] = f.refine_words.size();
                if (f.refine_dirty[j])
                    for (auto [w, mask] : layer_words[j]) {
                        if (!charge(1)) {
                            return false;
                        }
                        ull active = f.pv[w] & mask;
                        if (active) f.refine_words.emplace_back(w, active);
                    }
                f.refine_offsets[j + 1] = f.refine_words.size();
            }
            bool rows_removed = false;
            each(f.pu, [&](uint u) {
                for (uint j : f.refine_layers)
                    if (f.refine_dirty[j]) {
                        if (!charge(1)) {
                            return;
                        }
                        uint degree = selected_sizes[j];
                        for (size_t i = f.refine_offsets[j]; i < f.refine_offsets[j + 1]; ++i) {
                            auto [w, mask] = f.refine_words[i];
                            if (!charge(1)) {
                                return;
                            }
                            degree += __builtin_popcountll(rows[u][w] & mask);
                            if (degree >= beta[j]) break;
                        }
                        if (degree < beta[j]) {
                            erase(f.pu, u);
                            rows_removed = changed = true;
                            return;
                        }
                    }
            });
            if (refused) return false;
            if (!rows_removed || !check_columns) break;
            for (uint j : f.refine_layers) {
                if (!charge(1)) {
                    return false;
                }
                f.refine_dirty[j] = 0;
            }
            bool columns_removed = false;
            each(f.pv, [&](uint v) {
                uint degree = selected_sizes[p];
                for (size_t w = 0; w < wu && degree < beta[p]; ++w) {
                    if (!charge(1)) {
                        return;
                    }
                    degree += __builtin_popcountll(columns[v][w] & f.pu[w]);
                }
                if (degree < beta[p]) {
                    discard_column(v);
                    columns_removed = true;
                }
            });
            if (refused) return false;
            if (!columns_removed) break;
        }
        return !changed || promising(f, np, ns);
    }


    ull support_upper(Frame& f, uint capacity) {
        f.allowed.resize(wu);
        for (size_t w = 0; w < wu; ++w) {
            if (!charge(1)) {
                return 0;
            }
            f.allowed[w] = selected_u_bits[w] | f.pu[w];
        }
        const size_t stride = static_cast<size_t>(capacity) + 1;
        f.hist.resize(static_cast<size_t>(g.K) * stride);
        for (uint& cell : f.hist) {
            if (!charge(1)) {
                return 0;
            }
            cell = 0;
        }
        f.support.clear();
        f.support.reserve(f.right.size());
        for (size_t i = 0; i < f.right.size(); ++i) {
            uint v = f.right[i], degree;
            degree = intersection_count(columns[v], f.allowed);
            ++f.hist[secondary[v].first * stride + degree];
            f.support.push_back(degree);
        }
        f.sizes.assign(g.K, 0);
        if (!charge(g.K)) {
            return 0;
        }
        ull bound = 0;
        uint minimum = max(beta[p], static_cast<uint>(selected_u.size()));
        for (uint s = capacity;; --s) {
            bool feasible = true;
            f.sizes[p] = s;
            for (uint j = 0; j < g.K; ++j)
                if (j != p) {
                    if (!charge(1)) {
                        return 0;
                    }
                    f.sizes[j] += f.hist[j * stride + s];
                    uint cap = f.sizes[j];
                    if (cap < beta[j]) feasible = false;
                }
            if (feasible) {
                ull value = compute_star_size(f.sizes, objective, p);
                bound = max(bound, value);
            }
            if (s == minimum) break;
        }
        f.refine_support_ready = true;
        return bound;
    }

    bool promising(Frame& f, uint& np, uint& ns) {
        f.refine_support_ready = false;
        np = count(f.pu);
        ns = count(f.pv);
        f.sizes = selected_sizes;
        if (!charge(g.K)) {
            return false;
        }
        f.sizes[p] += np;
        f.right.clear();
        f.right.reserve(selected_v.size() + ns);
        for (uint v : selected_v) {
            if (!charge(1)) {
                return false;
            }
            f.right.push_back(v);
        }
        each(f.pv, [&](uint v) {
            ++f.sizes[secondary[v].first];
            f.right.push_back(v);
        });
        for (uint j = 0; j < g.K; ++j) {
            if (!charge(1)) {
                return false;
            }
            if (f.sizes[j] < beta[j]) return false;
        }
        ull upper = compute_star_size(f.sizes, objective, p);
        if ((!publish && upper <= best.size) || (improve && !improve(upper, secondary, f.right)))
            return false;
        if (refused) return false;
        if (!publish || improve) {
            upper = support_upper(f, f.sizes[p]);
            if (refused) return false;
            if ((!publish && upper <= best.size) || (improve && !improve(upper, secondary, f.right)))
                return false;
        }
        return true;
    }

    void emit(Frame& f, const Bits* add_u = nullptr, const Bits* add_v = nullptr) {
        f.sizes = selected_sizes;
        if (!charge(g.K)) {
            return;
        }
        if (add_u) f.sizes[p] += count(*add_u);
        f.right.clear();
        f.right.reserve(selected_v.size() + (add_v ? nv : 0));
        for (uint v : selected_v) {
            if (!charge(1)) {
                return;
            }
            f.right.push_back(v);
        }
        if (add_v)
            each(*add_v, [&](uint v) {
                ++f.sizes[secondary[v].first];
                f.right.push_back(v);
            });
        for (uint j = 0; j < g.K; ++j) {
            if (!charge(1)) {
                return;
            }
            if (f.sizes[j] < beta[j]) return;
        }

        ull value = compute_star_size(f.sizes, objective, p);
        bool offer = static_cast<bool>(publish);
        if (offer && improve) offer = improve(value, secondary, f.right);
        if (refused || (value <= best.size && !offer)) return;
        StarClique witness(g.K);
        witness.size = value;
        for (uint u : selected_u) {
            if (!charge(1)) {
                return;
            }
            witness.sets[p].push_back(primary[u]);
        }
        if (add_u) each(*add_u, [&](uint u) { witness.sets[p].push_back(primary[u]); });
        for (uint v : f.right) {
            if (!charge(1)) {
                return;
            }
            witness.sets[secondary[v].first].push_back(secondary[v].second);
        }
        if (refused) return;
        if (offer) {
            publish(witness);
            if (!charge(0)) {
                return;
            }
        }
        if (value > best.size) best = std::move(witness);
    }


    bool excluded_covers(const Bits& excluded, const Bits& added, const vector<Bits>& adjacency) {
        for (size_t w = 0; w < excluded.size(); ++w) {
            if (!charge(1)) {
                return false;
            }
            ull word = excluded[w];
            while (word) {
                uint v = static_cast<uint>(64 * w + __builtin_ctzll(word));
                word &= word - 1;
                bool covers = true;
                for (size_t k = 0; k < added.size(); ++k) {
                    if (!charge(1)) {
                        return false;
                    }
                    if (added[k] & ~adjacency[v][k]) {
                        covers = false;
                        break;
                    }
                }
                if (covers) return true;
            }
        }
        return false;
    }



    bool empty_cross_terminal(Frame& f, uint np, uint ns) {
        for (size_t w = 0; w < f.pu.size(); ++w) {
            if (!charge(1)) {
                return false;
            }
            ull word = f.pu[w];
            while (word) {
                uint u = static_cast<uint>(64 * w + __builtin_ctzll(word));
                word &= word - 1;
                for (size_t k = 0; k < wv; ++k) {
                    if (!charge(1)) {
                        return false;
                    }
                    if (rows[u][k] & f.pv[k]) return false;
                }
            }
        }

        if (!charge(0)) {
            return false;
        }

        if ((np || !ns) && empty(f.xu) && !excluded_covers(f.xv, f.pu, columns)) emit(f, &f.pu);
        if (ns && empty(f.xv) && !excluded_covers(f.xu, f.pv, rows)) emit(f, nullptr, &f.pv);
        return true;
    }



    void visit(size_t depth) {
        if (refused) return;

        ++nodes;
        if (!charge(1)) {
            return;
        }

        visit_body(depth);
    }



    // Branch on pivot non-neighbors to preserve complete enumeration.
    bool make_branches(Frame& f, uint np, uint ns) {

        uint pivot = 0, best_branches = UINT_MAX;
        bool primary_pivot = true;
        auto consider = [&](const Bits& set, bool left, bool in_p) {
            each(set, [&](uint v) {
                uint neighbors =
                    left ? intersection_count(rows[v], f.pv) : intersection_count(columns[v], f.pu);
                uint branches = (left ? ns : np) - neighbors + (in_p ? 1U : 0U);
                if (branches < best_branches) {
                    best_branches = branches;
                    pivot = v;
                    primary_pivot = left;
                }
            });
        };
        consider(f.xu, true, false);
        if (refused || !best_branches) return false;
        consider(f.xv, false, false);
        if (refused || !best_branches) return false;
        consider(f.pu, true, true);
        consider(f.pv, false, true);
        if (refused) return false;
        f.branch_u.clear();
        f.branch_v.clear();
        if (primary_pivot) {
            if ((f.pu[pivot / 64] >> (pivot % 64)) & 1ULL) f.branch_u.push_back(pivot);
            each(f.pv, [&](uint v) {
                if (!((rows[pivot][v / 64] >> (v % 64)) & 1ULL)) f.branch_v.push_back(v);
            });
        } else {
            each(f.pu, [&](uint u) {
                if (!((columns[pivot][u / 64] >> (u % 64)) & 1ULL)) f.branch_u.push_back(u);
            });
            if ((f.pv[pivot / 64] >> (pivot % 64)) & 1ULL) f.branch_v.push_back(pivot);
        }
        return true;
    }

    void visit_body(size_t depth) {
        maximum_primary = max(maximum_primary, selected_u.size());
        Frame& f = frames[depth];
        uint np = 0, ns = 0;
        if (!promising(f, np, ns) || refused) return;
        if (!refine_compact(f, np, ns) || refused) return;
        if (!np && !ns) {
            if (empty(f.xu) && empty(f.xv)) emit(f);
            return;
        }
        if (empty_cross_terminal(f, np, ns) || refused) return;

        if (!make_branches(f, np, ns) || refused) return;
        if (frames.size() <= depth + 1) frames.emplace_back();

        for (uint u : f.branch_u) {
            if (refused) return;
            Frame& child = frames[depth + 1];
            copy(child.pu, f.pu);
            erase(child.pu, u);
            copy(child.xu, f.xu);
            intersect(child.pv, f.pv, rows[u]);
            intersect(child.xv, f.xv, rows[u]);
            if (refused) return;
            selected_u.push_back(u);
            put(selected_u_bits, u);
            ++selected_sizes[p];
            visit(depth + 1);
            --selected_sizes[p];
            erase(selected_u_bits, u);
            selected_u.pop_back();
            erase(f.pu, u);
            put(f.xu, u);
        }
        for (uint v : f.branch_v) {
            if (refused) return;
            Frame& child = frames[depth + 1];
            intersect(child.pu, f.pu, columns[v]);
            intersect(child.xu, f.xu, columns[v]);
            copy(child.pv, f.pv);
            erase(child.pv, v);
            copy(child.xv, f.xv);
            if (refused) return;
            selected_v.push_back(v);
            ++selected_sizes[secondary[v].first];
            visit(depth + 1);
            --selected_sizes[secondary[v].first];
            selected_v.pop_back();
            erase(f.pv, v);
            put(f.xv, v);
        }
    }

    void prepare() {
        if (prepared) return;
        prepared = true;
        selected_u.reserve(nu);
        selected_v.reserve(nv);
        if (!charge(nu + nv + g.K + wu)) {
            return;
        }

        {
            layer_words.resize(g.K);
            if (!charge(g.K)) {
                return;
            }
            for (size_t v = 0; v < nv; ++v) {
                if (!charge(1)) {
                    return;
                }
                auto& words = layer_words[secondary[v].first];
                if (words.empty() || words.back().first != v / 64) words.emplace_back(v / 64, 0);
                words.back().second |= 1ULL << (v % 64);
            }
        }
    }

  public:
    ExactPivotEnumeration(
        const KPartiteGraph& graph, uint primary_layer, ProblemType problem, const vector<uint>& thresholds,
        const vector<uint>& primary_vertices, const vector<pair<uint, uint>>& secondary_vertices,
        const vector<Bits>& bit_rows, const vector<Bits>& bit_columns, StarClique& best_result,
        const function<void(const StarClique&)>& callback,
        const function<bool(ull, const vector<pair<uint, uint>>&, const vector<uint>&)>& can_improve,
        ull& visited, Charge budget_charge, bool& budget_refused)
        : g(graph), p(primary_layer), objective(problem), beta(thresholds), primary(primary_vertices),
          secondary(secondary_vertices), rows(bit_rows), columns(bit_columns), best(best_result),
          publish(callback), improve(can_improve), nodes(visited), charge(std::move(budget_charge)),
          refused(budget_refused), nu(primary.size()), nv(secondary.size()), wu((nu + 63) / 64),
          wv((nv + 63) / 64), selected_sizes(g.K, 0), selected_u_bits(wu, 0) {}



    void run() {
        frames.emplace_back();
        initialize(frames[0]);
        if (refused) return;
        prepare();
        if (!refused) visit(0);
    }
};

enum class BudgetStatus { Complete, Partial, MemoryFallback };

class ExactCenterClosure {
    const KPartiteGraph& g;
    uint p, q, a, b;
    ProblemType objective;
    const vector<uint>& beta;
    vector<uint> primary;
    vector<pair<uint, uint>> secondary;
    vector<vector<ull>> rows, columns;
    StarClique best;
    bool all_graph;
    function<void(const StarClique&)> publish;
    const vector<uint>* primary_domain;
    uint anchor;
    function<bool(ull, const vector<pair<uint, uint>>&, const vector<uint>&)> improve;

    ull nodes = 0;
    ull work = 0;
    ull next_budget_check = 0, last_checked_node = ULLONG_MAX;
    function<bool(ull, ull)> keep_going;
    bool refused = false;

    size_t nw = 0, nr = 0;

    // Account work and poll the budget without discarding unresolved coverage.
    template <bool limited> bool charge(ull amount) {
        if (refused) return false;
        if constexpr (limited) {
            work = amount > ULLONG_MAX - work ? ULLONG_MAX : work + amount;

            if (!amount || work >= next_budget_check || nodes != last_checked_node) {
                last_checked_node = nodes;
                next_budget_check = work > ULLONG_MAX - 64 ? ULLONG_MAX : work + 64;
                if (keep_going && !keep_going(nodes, work)) refused = true;
            }
        }
        return !refused;
    }

    static bool has(const vector<ull>& bits, size_t x) { return (bits[x / 64] >> (x % 64)) & 1ULL; }

    static void put(vector<ull>& bits, size_t x) { bits[x / 64] |= 1ULL << (x % 64); }

    static uint count(const vector<ull>& bits) {
        uint n = 0;
        for (ull x : bits)
            n += __builtin_popcountll(x);
        return n;
    }

  public:
    ull visited_nodes() const { return nodes; }


    ExactCenterClosure(
        const KPartiteGraph& graph, ProblemType pt, const vector<uint>& thresholds, uint pp, uint u,
        uint secondary_part, uint v, bool full = false, function<void(const StarClique&)> callback = {},
        const vector<uint>* domain = nullptr, uint primary_anchor = UINT_MAX,
        function<bool(ull, const vector<pair<uint, uint>>&, const vector<uint>&)> can_improve = {})
        : g(graph), p(pp), q(secondary_part), a(u), b(v), objective(pt), beta(thresholds), best(graph.K),
          all_graph(full), publish(std::move(callback)), primary_domain(domain), anchor(primary_anchor),
          improve(std::move(can_improve)) {}


    BudgetStatus enumerate_budgeted(StarClique& result, function<bool(ull, ull)> budget) {
        keep_going = std::move(budget);
        bool complete = solve_impl<true>(result);
        keep_going = {};
        if (refused) return BudgetStatus::Partial;
        return complete ? BudgetStatus::Complete : BudgetStatus::MemoryFallback;
    }

  private:
    // Prepare typed adjacency matrices, then enumerate within the admitted budget.
    template <bool limited> bool solve_impl(StarClique& result) {
        if (primary_domain) {
            if (!charge<limited>(primary_domain->size())) return false;
            primary = *primary_domain;
            for (uint j = 0; j < g.K; ++j)
                if (j != p)
                    for (uint v : g.neighbors(p, anchor, j))
                        secondary.emplace_back(j, v);
        } else if (all_graph) {
            if (!charge<limited>(g.vertex_ids[p].size())) return false;
            for (uint u = 0; u < g.vertex_ids[p].size(); ++u)
                primary.push_back(u);
            for (uint j = 0; j < g.K; ++j)
                if (j != p)
                    for (uint v = 0; v < g.vertex_ids[j].size(); ++v)
                        secondary.emplace_back(j, v);
        } else {
            if (!charge<limited>(g.neighbors(q, b, p).size())) return false;
            primary = g.neighbors(q, b, p);
            for (uint j = 0; j < g.K; ++j)
                if (j != p)
                    for (uint v : g.neighbors(p, a, j))
                        secondary.emplace_back(j, v);
        }
        nw = (primary.size() + 63) / 64;
        nr = (secondary.size() + 63) / 64;
        if (!charge<limited>(secondary.size())) return false;
        const size_t words_budget = (64ULL << 20) / sizeof(ull);
        auto fits = [&] {
            if (nw && secondary.size() > words_budget / nw) return false;
            size_t used = nw * secondary.size();
            return !nr || primary.size() <= (words_budget - used) / nr;
        };
        if (!fits()) {

            if (!fits()) return false;
        }
        vector<uint> sizes(g.K, 0);
        sizes[p] = static_cast<uint>(primary.size());
        for (auto [j, v] : secondary) {
            (void)v;
            ++sizes[j];
        }
        for (uint j = 0; j < g.K; ++j)
            if (sizes[j] < beta[j]) {
                result = best;
                return true;
            }
        if (!charge<limited>(primary.size() * nr + secondary.size() * nw)) return false;
        rows.assign(primary.size(), vector<ull>(nr, 0));
        columns.assign(secondary.size(), vector<ull>(nw, 0));

        for (size_t x = 0; x < primary.size(); ++x) {
            size_t begin = 0;
            for (uint j = 0; j < g.K; ++j)
                if (j != p) {
                    size_t end = begin;
                    while (end < secondary.size() && secondary[end].first == j)
                        ++end;
                    const auto& neighbors = g.neighbors(p, primary[x], j);
                    if (!charge<limited>(end - begin + neighbors.size())) return false;
                    size_t i = begin, k = 0;
                    while (i < end && k < neighbors.size()) {
                        uint v = secondary[i].second;
                        if (v < neighbors[k])
                            ++i;
                        else if (v > neighbors[k])
                            ++k;
                        else {
                            put(rows[x], i);
                            put(columns[i], x);
                            ++i;
                            ++k;
                        }
                    }
                    begin = end;
                }
        }
        {
            if (!charge<limited>(primary.size() * (nr + secondary.size() + 16))) return false;
            vector<size_t> order(primary.size());
            iota(order.begin(), order.end(), 0);
            vector<uint> degree(primary.size());
            for (size_t x = 0; x < primary.size(); ++x)
                degree[x] = count(rows[x]);
            stable_sort(order.begin(), order.end(),
                        [&](size_t x, size_t y) { return degree[x] < degree[y]; });
            auto old_primary = primary;
            auto old_rows = std::move(rows);
            rows.resize(primary.size());
            for (auto& col : columns)
                fill(col.begin(), col.end(), 0);
            for (size_t x = 0; x < primary.size(); ++x) {
                primary[x] = old_primary[order[x]];
                rows[x] = std::move(old_rows[order[x]]);
                for (size_t w = 0; w < rows[x].size(); ++w) {
                    ull bits = rows[x][w];
                    while (bits) {
                        size_t v = 64 * w + __builtin_ctzll(bits);
                        bits &= bits - 1;
                        put(columns[v], x);
                    }
                }
            }
        }

        auto account = [&](ull amount) { return charge<limited>(amount); };
        ExactPivotEnumeration<decltype(account)> enumeration(g, p, objective, beta, primary, secondary, rows,
                                                             columns, best, publish, improve, nodes, account,
                                                             refused);
        enumeration.run();
        result = std::move(best);
        return !refused;
    }
};
} // namespace mcviz::legacy_star
