#include "pattern_partner.h"
#include "common/text.h"
#include <algorithm>
#include <limits>
#include <numeric>

namespace mcviz {
namespace {
using Word = std::uint64_t;
constexpr auto missing = std::numeric_limits<std::size_t>::max();
}
// Directed required-layer blocks use S=sum_(i,j) n_i*ceil(n_j/64) words.
// Construction O(N+M_P+S), storage O(N+S); budget refusal preserves all domains.
PartnerRows::PartnerRows(const Graph& g, const LayerPattern& p, const Sets& base, std::uint64_t memory)
    : g_(g), p_(p) {
    long double estimate = sizeof(*this) + 128.0L * g.K;
    for (Layer i = 0; i < g.K; ++i) {
        estimate += 8.0L * g.vertex_ids[i].size() + 4.0L * base[i].size();
        for (auto j : p.neighbors[i]) estimate += 64.0L + 8.0L * base[i].size() * ((base[j].size()+63)/64);
    }
    if (estimate > memory) return;
    {
        base_ = base;

        pos_.resize(g.K); blocks_.resize(g.K);
        for (Layer i = 0; i < g.K; ++i) {
            pos_[i].assign(g.vertex_ids[i].size(), missing);

            for (std::size_t a = 0; a < base[i].size(); ++a) pos_[i][base[i][a]] = a;
        }
        for (Layer i = 0; i < g.K; ++i) for (auto j : p.neighbors[i]) {
            Block b{j, (base[j].size()+63)/64, {}};
            b.bits.assign(base[i].size()*b.words, 0);
            for (std::size_t a = 0; a < base[i].size(); ++a) for (auto v : g.neighbors(i,base[i][a],j)) {
                const auto c = pos_[j][v];
                if (c != missing) b.bits[a*b.words+c/64] |= Word(1)<<(c%64);
            }

            blocks_[i].push_back(std::move(b));
        }
        bytes_ = static_cast<std::uint64_t>(estimate); ready_ = true;
    }
}

// In any feasible tuple containing u, its beta_i-1 same-layer partners each
// share at least beta_j neighbors with u in EVERY required neighbor layer j.
// Delete only after exhausting candidate partners without enough successes.
// An interrupted pair test retains u. Sequential batches plus queue-core reach a
// fixed point unless the work budget runs out; completed deletions remain safe.
// Per pass O(sum_i n_i^2 sum_j ceil(base_j/64)); O(sum_j ceil(base_j/64)+N)
// scratch plus queue core. No subset/pair table is stored, no exponential space.
CoreResult PartnerRows::reduce(const Thresholds& beta, Sets domains, const Required& fixed,
                              std::uint64_t budget, PartnerStats& stats, CoreWorkspace& workspace) const {
    CoreResult result{std::move(domains), true, 0, 0};
    if (!ready_) {
        result.possible = enough_capacity(result.domains,beta) && contains_required(result.domains,fixed);
        return result;
    }
    while (true) {
        ++stats.rounds;
        std::vector<std::vector<Word>> masks(g_.K);

        for (Layer i = 0; i < g_.K; ++i) {
            masks[i].assign((base_[i].size()+63)/64, 0);
            for (auto u : result.domains[i]) {
                const auto a = pos_[i][u];
                masks[i][a/64] |= Word(1)<<(a%64);
            }
        }

        Sets next = result.domains;

        bool stopped = false, changed = false;
        for (Layer i = 0; i < g_.K && !stopped; ++i) {
            if (beta[i] <= 1) continue;
            std::vector<std::size_t> order(blocks_[i].size());
            std::iota(order.begin(),order.end(),0);
            std::sort(order.begin(),order.end(),[&](auto a,auto b) {
                return result.domains[blocks_[i][a].target].size() < result.domains[blocks_[i][b].target].size();
            });

            next[i].clear();
            for (std::size_t uidx=0; uidx<result.domains[i].size(); ++uidx) {
                auto u = result.domains[i][uidx];
                std::uint64_t partners = 0;

                for (auto v : result.domains[i]) {
                    if (u==v) continue;
                    if (stats.work >= budget) { stopped=true; break; }

                    ++stats.work;
                    bool compatible = true;

                    for (auto bidx : order) {
                        const auto& b = blocks_[i][bidx];
                        const auto au=pos_[i][u]*b.words, av=pos_[i][v]*b.words;
                        std::uint64_t common=0;

                        for (std::size_t z=0; z<b.words && common<beta[b.target]; ++z) {
                            if (stats.work >= budget) { stopped=true; break; }

                            ++stats.work;

                            common += __builtin_popcountll(b.bits[au+z]&b.bits[av+z]&masks[b.target][z]);
                        }
                        if (stopped || common<beta[b.target]) { compatible=false; break; }
                    }
                    if (stopped) break;
                    if (compatible && ++partners >= beta[i]-1) break;
                }
                if (stopped) {
                    next[i].insert(next[i].end(),result.domains[i].begin()+uidx,result.domains[i].end());

                    break;
                }

                if (partners >= beta[i]-1) next[i].push_back(u);
                else { ++stats.removed; ++result.removed; changed=true; }
            }
        }

        if (stopped) ++stats.budget_stops;

        result.domains = std::move(next);

        if (changed) {
            auto core = mcviz::reduce(g_,p_,beta,std::move(result.domains),fixed,CoreBackend::queue,&workspace);
            stats.core_removed += core.removed;
            result.removed += core.removed;
            result.domains = std::move(core.domains);
            result.possible = core.possible;
        }
        if (!result.possible || !changed || stopped) break;
    }

    result.possible = result.possible && enough_capacity(result.domains,beta) && contains_required(result.domains,fixed);
    return result;
}
} // namespace mcviz
