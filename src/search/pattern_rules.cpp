#include "pattern_rules.h"
#include "common/text.h"
#include <algorithm>
#include <optional>

namespace mcviz {
namespace {
bool spend(std::uint64_t& left, std::uint64_t a, std::uint64_t b = 1) {
    if (a && b > left / a) return false;
    left -= a * b;
    return true;
}

bool estimate(const LayerPattern& p, const Sets& d, std::uint64_t& left) {
    if (!spend(left, p.K)) return false;
    for (auto e : p.edges) {
        if (!spend(left, d[e.first].size(), d[e.second].size())) return false;
        if (!spend(left, d[e.first].size(), d[e.second].size())) return false;
        if (!spend(left, d[e.first].size(), d[e.second].size())) return false;
        if (!spend(left, d[e.first].size(), d[e.second].size())) return false;
        for (auto i : {e.first, e.second}) {
            std::uint64_t bits = 1;
            for (auto n = d[i].size(); n > 1; n >>= 1) ++bits;
            if (!spend(left, d[i].size(), bits)) return false;
        }
    }
    return true;
}
struct PairDegrees {
    std::vector<Score> parent, child;
    Score mandatory_parent = 0, mandatory_child = 0;
};

std::optional<Score> tree_upper(const Graph& g, const LayerPattern& p, Objective obj,
                                const Thresholds& lower, const Sets& d, const Sets& q,
                                const std::vector<Layer>& parent, const std::vector<Layer>& order) {
    std::vector<PairDegrees> pairs(p.K);
    for (Layer v = 1; v < p.K; ++v) {
        Layer u = parent[v];
        auto& a = pairs[v];
        a.parent.assign(d[u].size(), 0); a.child.assign(d[v].size(), 0);
        a.mandatory_parent = d[v].size(); a.mandatory_child = d[u].size();
        for (std::size_t i = 0; i < d[u].size(); ++i) for (std::size_t j = 0; j < d[v].size(); ++j) {
            if (g.has_edge(u, d[u][i], v, d[v][j])) { ++a.parent[i]; ++a.child[j]; }
        }
        for (std::size_t i = 0; i < d[u].size(); ++i) if (std::binary_search(q[u].begin(), q[u].end(), d[u][i]))
            a.mandatory_parent = std::min(a.mandatory_parent, a.parent[i]);
        for (std::size_t i = 0; i < d[v].size(); ++i) if (std::binary_search(q[v].begin(), q[v].end(), d[v][i]))
            a.mandatory_child = std::min(a.mandatory_child, a.child[i]);
        std::sort(a.parent.rbegin(), a.parent.rend()); std::sort(a.child.rbegin(), a.child.rend());
    }
    std::vector<std::vector<std::optional<Score>>> dp(p.K);
    for (std::size_t at = order.size(); at-- > 0;) {
        Layer u = order[at];
        dp[u].resize(d[u].size() + 1);
        for (std::size_t a = static_cast<std::size_t>(lower[u]); a <= d[u].size(); ++a) {
            std::optional<Score> value = obj == Objective::edge ? 0 : a;
            for (auto v : p.neighbors[u]) if (parent[v] == u) {
                std::optional<Score> best;
                const auto& deg = pairs[v];
                for (std::size_t b = static_cast<std::size_t>(lower[v]); b <= d[v].size(); ++b) {
                    if (!dp[v][b] || deg.parent[a - 1] < b || deg.child[b - 1] < a ||
                        deg.mandatory_parent < b || deg.mandatory_child < a) continue;
                    Score candidate = *dp[v][b];
                    if (obj == Objective::edge) candidate = checked_add(candidate, checked_mul(a, b));
                    if (!best || candidate > *best) best = candidate;
                }
                if (!best) { value.reset(); break; }
                if (obj == Objective::balanced) value = std::min(*value, *best);
                else value = checked_add(*value, *best);
            }
            dp[u][a] = value;
        }
    }
    std::optional<Score> best;
    for (auto v : dp[0]) if (v && (!best || *v > *best)) best = v;
    return best;
}
}

RulesUpperResult reviewed_upper(const Graph& g, const LayerPattern& p, Objective obj,
                                const Thresholds& lower_input, const Sets& domains, const Required& required,
                                std::uint64_t& budget, CoreBackend core_backend,
                                CoreWorkspace* core_workspace) {
    auto conditioned = condition_domains(g, p, domains, required);
    const auto& d = conditioned.domains;
    Thresholds lower = lower_input;
    Sets q(p.K);
    for (auto v : required) q[v.layer].push_back(v.pos);
    for (Layer i = 0; i < p.K; ++i) {
        std::sort(q[i].begin(), q[i].end()); q[i].erase(std::unique(q[i].begin(), q[i].end()), q[i].end());
        lower[i] = std::max<Score>(lower[i], q[i].size());
    }
    RulesUpperResult out;
    out.possible = conditioned.possible && enough_capacity(d, lower);
    if (!out.possible) return out;
    out.upper = capacity_bound(p, obj, d);
    if (p.edges.empty()) return out;
    std::vector<Layer> parent(p.K, p.K), order{0};
    parent[0] = 0;
    for (std::size_t j = 0; j < order.size(); ++j) for (auto v : p.neighbors[order[j]]) if (parent[v] == p.K) {
        parent[v] = order[j]; order.push_back(v);
    }
    if (order.size() == p.K && p.edges.size() == p.K - 1) {
        auto left = budget;
        if (estimate(p, d, left)) {
            budget = left; out.tree_used = true;
            auto value = tree_upper(g, p, obj, lower, d, q, parent, order);
            if (!value) { out.possible = false; out.upper = 0; return out; }
            out.upper = std::min(out.upper, *value);
        } else out.budget_fallback = true;
    }
    if (obj == Objective::balanced) {
        Score cap = capacity_bound(p, Objective::balanced, d);
        std::uint64_t n = 0, trials = 2;
        for (const auto& a : d) n += a.size();
        for (auto x = cap; x > 1; x >>= 1) ++trials;
        auto left = budget;
        bool fits = true;

        auto one_left = budget;
        if (!estimate(p, d, one_left)) fits = false;
        std::uint64_t one = budget - one_left;
        if (fits && (!spend(left, n + 1, trials) || !spend(left, one, (n + 1) * trials))) fits = false;
        if (!fits) out.budget_fallback = true;
        else {
            budget = left; out.core_used = true;
            if (!reduce(g, p, lower, d, required, core_backend, core_workspace).possible) {
                out.possible = false; out.upper = 0; return out;
            }
            Score lo = 0, hi = cap;
            while (lo < hi) {
                Score mid = lo + (hi - lo) / 2 + (hi - lo) % 2;
                auto t = lower;
                for (auto& v : t) v = std::max(v, mid);
                if (reduce(g, p, t, d, required, core_backend, core_workspace).possible) lo = mid;
                else hi = mid - 1;
            }
            out.upper = std::min(out.upper, lo);
        }
    }
    return out;
}

Thresholds improving_thresholds(const LayerPattern& p, Objective obj, const Sets& d,
                                const Thresholds& gamma, const Sets& mandatory, Score L) {
    Thresholds lower = gamma;
    Score U = capacity_bound(p, obj, d);

    for (Layer i = 0; i < p.K; ++i) {
        Score rho = 0;
        if (obj == Objective::edge) {
            Score c = 0;
            for (auto j : p.neighbors[i]) c = checked_add(c, d[j].size());
            Score a = U - checked_mul(d[i].size(), c);
            if (c && L >= a) rho = checked_add((L - a) / c, 1);
        } else if (obj == Objective::vertex) {
            Score other = U - d[i].size();
            if (L >= other) rho = checked_add(L - other, 1);
        } else rho = checked_add(L, 1);
        lower[i] = std::max({lower[i], Score(mandatory[i].size()), rho});
    }
    return lower;
}
}
