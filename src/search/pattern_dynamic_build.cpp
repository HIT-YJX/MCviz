#include "pattern_dynamic_build.h"
#include "pattern_adjacency.h"
#include "pattern_edge_domains.h"
#include "pattern_partner.h"
#include "pattern_rules.h"
#include "common/text.h"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <ctime>
#include <map>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <tuple>

namespace mcviz {
namespace {
using Clock=std::chrono::steady_clock;
double ms(Clock::time_point start) { return std::chrono::duration<double,std::milli>(Clock::now()-start).count(); }
std::uint64_t vertices(const Sets& d) { std::uint64_t n=0; for (const auto& x:d) n+=x.size(); return n; }
std::uint64_t key_bytes(const Sets& d) {
    std::uint64_t n=128+sizeof(Sets)+d.capacity()*sizeof(std::vector<VertexPos>);
    for (const auto& x:d) n+=x.capacity()*sizeof(VertexPos);
    return n;
}
struct Worker {
    std::atomic<std::uint64_t> completed{0};
    std::uint64_t calls=0,nodes=0,zero=0,hits=0,entries=0,cache_bytes=0,cache_refused=0;
    std::uint64_t local_core_calls=0,local_core_removed=0,local_core_scans=0,sparse_requests=0;
    std::uint64_t certified=0,warm=0,warm_positive=0,positive=0,no_improvement=0,published=0,publications=0;
    std::uint64_t improvement_core=0,improvement_raised=0,improvement_rejected=0;
    std::uint64_t cache_lookups=0,cache_duplicates=0,cache_foreign_hits=0;
    PartnerStats partner;
    SolverStats solver_stats;
    double edge_ms=0,solver_ms=0,partner_ms=0,domain_ms=0,publish_ms=0,improvement_ms=0,key_ms=0;
};

// Hash only routes a complete typed domain to one of 64 locks; ordered map
// equality still compares every layer/position. O(|D|+|D| log(Q_s+1)) per
// lookup/insert, where Q_s is this shard's entry count, plus lock contention.
// No solver or witness publication executes while a shard lock is held.
class SharedScores {
    struct Entry { Score value; unsigned owner; };
    struct Shard { std::mutex mutex; std::map<Sets,Entry> values; };
    static constexpr std::size_t count=64;
    std::unique_ptr<Shard[]> shards_;
    std::uint64_t budget_=0;
    std::atomic<std::uint64_t> used_{0};
public:
    enum class Insert { added, duplicate, refused };
    explicit SharedScores(std::uint64_t budget):budget_(budget) {
        const auto overhead=sizeof(*this)+count*sizeof(Shard);
        if(budget<overhead) return;
        { shards_=std::make_unique<Shard[]>(count);used_.store(overhead); }
    }
    bool ready() const { return bool(shards_); }
    std::uint64_t bytes() const { return used_.load(std::memory_order_relaxed); }
    std::size_t size() const { return ready()?count:0; }
    // Layer sizes delimit the original-position lists. Collisions affect only
    // lock sharing, never equivalence; no key storage or allocation is needed.
    std::size_t shard(const Sets& d) const {
        std::uint64_t h=1469598103934665603ULL;
        auto mix=[&](std::uint64_t x) { h^=x;h*=1099511628211ULL; };
        mix(d.size());
        for(const auto& row:d) { mix(row.size());for(auto u:row) mix(u); }
        h^=h>>33;h*=0xff51afd7ed558ccdULL;h^=h>>33;
        return h%count;
    }
    bool get(std::size_t s,const Sets& d,Score& value,unsigned owner,bool& foreign) {
        if(!ready()) return false;
        auto& part=shards_[s];
        std::lock_guard<std::mutex> lock(part.mutex);
        const auto at=part.values.find(d);
        if(at==part.values.end()) return false;
        value=at->second.value;foreign=at->second.owner!=owner;
        return true;
    }
    // Reserve from one total logical byte budget while holding one shard lock.
    // Duplicate completion consumes no cache budget.
    // No waiting dependency or partially solved cache entry is ever installed.
    Insert put(std::size_t s,Sets&& d,Score value,unsigned owner) {
        if(!ready()) return Insert::refused;
        auto& part=shards_[s];
        std::lock_guard<std::mutex> lock(part.mutex);
        const auto at=part.values.find(d);
        if(at!=part.values.end()) {
            return Insert::duplicate;
        }
        const auto bytes=key_bytes(d);
        auto used=used_.load(std::memory_order_relaxed);
        do {
            if(bytes>budget_-used) return Insert::refused;
        } while(!used_.compare_exchange_weak(used,used+bytes,std::memory_order_relaxed));
        { part.values.emplace(std::move(d),Entry{value,owner}); }
        return Insert::added;
    }
};

// Publish only actual witness scores to their contained original P-edges.
// Relaxed atomic max transfers a number, never a witness object or completion
// flag. O(|E_P(W)| log(M+1)) edge lookup plus contention retries, O(1) scratch.
void publish(const Graph& g,const LayerPattern& p,const SearchResult& result,
             std::vector<std::atomic<Score>>& lower,Worker& w) {
    if(!result.found) return;
    const auto start=Clock::now();
    auto less=[](const GraphEdge& a,const GraphEdge& b) {
        return std::tie(a.u.layer,a.u.pos,a.v.layer,a.v.pos)<std::tie(b.u.layer,b.u.pos,b.v.layer,b.v.pos);
    };
    for(auto ij:p.edges) for(auto u:result.sets[ij.first]) for(auto v:result.sets[ij.second]) {
        GraphEdge e{{ij.first,u},{ij.second,v}};
        const auto at=std::lower_bound(g.edges.begin(),g.edges.end(),e,less);
        auto& bound=lower[at-g.edges.begin()];
        auto old=bound.load(std::memory_order_relaxed);
        while(old<result.score && !bound.compare_exchange_weak(old,result.score,std::memory_order_relaxed)) {}
        ++w.publications;
    }
    ++w.published;w.publish_ms+=ms(start);
}
}

// Dynamic fetch-add assigns one original P-edge at a time, so no worker waits
// behind another worker's fixed slice. Shared adjacency is immutable; all domains,
// mandatory endpoints, solver stacks and witnesses are worker-local. Fused
// builds share completed exact scores.
// O(global reductions + max dynamic worker load + M) time excluding forest; exact
// search remains exponential. Space: shared O(N+bit blocks+M), plus T bounded
// caches and private core/solver workspaces. Cache budget is never multiplied by T.
PatternReuseBuildResult build_pattern_dynamic(const Graph& g,const LayerPattern& p,
    const Thresholds& beta,Sets& base,unsigned threads,std::uint64_t memory,std::uint64_t partner_work,
    std::uint64_t pivot_memory_bytes) {
    const auto started=Clock::now();
    const auto cpu_started=std::clock();
    PatternReuseBuildResult out;
    out.value_kind="exact"; out.values.resize(g.edges.size());
    out.counters["dynamic_fusion_enabled"]=true;
    out.counters["dynamic_beta_reduction_enabled"]=true;
    out.counters["dynamic_partner_enabled"]=true;
    out.counters["dynamic_adjacency_enabled"]=true;
    out.counters["dynamic_score_cache_enabled"]=true;
    out.counters["pattern_warm_cutoff_enabled"]=true;
    out.counters["pattern_certificate_enabled"]=true;
    out.counters["pattern_improvement_core_enabled"]=true;
    const auto score_budget=memory-memory/4-memory/4;
    SharedScores scores(score_budget);
    out.counters["dynamic_shared_cache_enabled"]=true;
    out.counters["dynamic_shared_cache_shards"]=scores.size();
    out.counters["dynamic_shared_cache_budget_bytes"]=score_budget;
    std::vector<std::atomic<Score>> lower(g.edges.size());
    for(auto& x:lower) x.store(0,std::memory_order_relaxed);
    auto begin=Clock::now();
    // Two immutable bit representations and all mutable caches share the budget.
    // Each representation receives one quarter; remaining half is shared for
    // completed exact scores.
    PartnerRows partners(g,p,base,memory/4);
    CoreWorkspace global_workspace;
    PartnerStats global;
    CoreResult reduced;
    reduced.possible=enough_capacity(base,beta);
    {
        reduced=partners.reduce(beta,base,{},partner_work,global,global_workspace);
        base=std::move(reduced.domains);
    }
    out.counters["dynamic_partner_ready"]=partners.ready();
    out.counters["dynamic_global_partner_work"]=global.work;
    out.counters["dynamic_global_partner_removed"]=global.removed;
    out.counters["dynamic_global_partner_core_removed"]=global.core_removed;
    out.counters["dynamic_global_partner_budget_stops"]=global.budget_stops;
    out.counters["dynamic_global_vertices"]=vertices(base);
    out.timing["dynamic_global_partner_ms"]=ms(begin);

    auto shared=std::make_shared<const PatternAdjacency>(g,p,base,true,memory/4);
    out.counters["dynamic_adjacency_ready"]=shared->ready();
    out.counters["dynamic_shared_bytes"]=partners.bytes()+shared->logical_bytes();
    out.timing["dynamic_adjacency_prepare_ms"]=shared->prepare_ms();
    std::vector<std::vector<unsigned char>> in(g.K);
    for (Layer i=0;i<g.K;++i) { in[i].resize(g.vertex_ids[i].size()); for(auto u:base[i]) in[i][u]=1; }
    std::vector<std::size_t> ids;
    for (std::size_t id=0;id<g.edges.size();++id) {
        const auto& e=g.edges[id];
        if(!p.adjacent(e.u.layer,e.v.layer)) continue;
        ++out.counters["pattern_edges"];
        if(reduced.possible && in[e.u.layer][e.u.pos] && in[e.v.layer][e.v.pos]) ids.push_back(id);
    }
    out.counters["dynamic_tasks"]=ids.size();
    out.counters["pattern_global_core_skipped"]=out.counters["pattern_edges"]-ids.size();
    out.timing["dynamic_prepare_ms"]=ms(started);
    const auto compute_started=Clock::now();
    std::vector<Worker> workers(threads);
    std::atomic<std::size_t> next{0};
    auto run=[&](unsigned t) {
        auto& w=workers[t];
        {
            PatternEdgeDomains domains(g,p,beta,base,true,true,memory/4,shared);
            CoreWorkspace workspace;
            for(;;) {
                const auto at=next.fetch_add(1);
                if(at>=ids.size()) break;
                const auto id=ids[at]; const auto& edge=g.edges[id];
                const Required fixed{edge.u,edge.v};
                 
                const auto edge_start=Clock::now();
                
                auto phase=Clock::now();
                auto d=domains.get(edge); w.domain_ms+=ms(phase);
                Score value=0; bool hit=false;
                double solver_ms=0;
                std::uint64_t nodes=0;
                if(d.possible) {
                     phase=Clock::now(); PartnerStats ps;
                    // Work limits apply only to optional necessary-condition filtering.
                    // The exact solver itself has no time/node budget or silent fallback.
                    d=partners.reduce(beta,std::move(d.domains),fixed,std::min<std::uint64_t>(partner_work,200000),ps,workspace);
                    w.partner_ms+=ms(phase); w.partner.work+=ps.work; w.partner.removed+=ps.removed;
                    w.partner.core_removed+=ps.core_removed; w.partner.budget_stops+=ps.budget_stops;
                    w.partner.rounds+=ps.rounds;
                }
                if(d.possible) {
                    const Score incumbent=lower[id].load(std::memory_order_relaxed);
                    w.positive+=incumbent>0;
                    // D remains endpoint-conditioned and includes both endpoints.
                    // Any feasible tuple in D can add them, so optimum(D) equals
                    // the anchored optimum. Exact full-D keys safely share values.
                    std::size_t shard=0;
                    bool foreign=false;
                    {
                        const auto key_start=Clock::now();
                        shard=scores.ready()?scores.shard(d.domains):0;
                        ++w.cache_lookups;
                        hit=scores.get(shard,d.domains,value,t,foreign);
                        w.key_ms+=ms(key_start);
                    }
                    if(hit) {
                        ++w.hits;w.cache_foreign_hits+=foreign;
                    }
                    else {
                        const Score upper=capacity_bound(p,Objective::edge,d.domains);
                        bool resolved=incumbent==upper;
                        if(resolved) { value=incumbent;++w.certified; }
                        Sets improved;
                        Thresholds thresholds=beta;
                        bool restricted=false;
                        if(!resolved) {
                            phase=Clock::now();
                            Sets mandatory(g.K);
                            for(auto q:fixed) mandatory[q.layer].push_back(q.pos);
                            thresholds=improving_thresholds(p,Objective::edge,d.domains,beta,mandatory,incumbent);
                            restricted=thresholds!=beta;
                            if(restricted) {
                                ++w.improvement_raised;++w.improvement_core;
                                auto r=reduce_queue(g,p,thresholds,d.domains,fixed,&workspace,nullptr);
                                resolved=!r.possible || capacity_bound(p,Objective::edge,r.domains)<=incumbent;
                                improved=std::move(r.domains);
                                if(resolved) { value=incumbent;++w.improvement_rejected; }
                            }
                            w.improvement_ms+=ms(phase);
                        }
                        if(!resolved) {
                             phase=Clock::now();
                            
                            SolverOptions options; options.reviewed_rules=true; options.research_skip_initial_core=true;
                            options.backend=RegionBackend::pivot_ordered;options.pivot_memory_bytes=pivot_memory_bytes;
                            options.core_backend=CoreBackend::queue; options.core_workspace=&workspace;
                            const auto cutoff=std::optional<Score>(incumbent);
                            auto r=solve_region(g,p,Objective::edge,thresholds,restricted?improved:d.domains,fixed,cutoff,options);
                            value=std::max(incumbent,r.best.found?r.best.score:0); nodes=r.stats.nodes;
                            ++w.calls;w.nodes+=nodes;solver_ms=ms(phase);w.solver_ms+=solver_ms;
                            {
                                ++w.warm;w.warm_positive+=incumbent>0;w.no_improvement+=!r.best.found;
                                publish(g,p,r.best,lower,w);
                            }
                            
                            merge_solver_stats(w.solver_stats,r.stats);
                        }
                        // The ORIGINAL beta domain supports endpoint augmentation.
                        // L has a genuine edge witness preserved in this domain;
                        // complete search of all >L solutions proves max(L,best)
                        // is its unanchored optimum too. Never key the PA-09 domain
                        // or cache an arbitrary caller cutoff/no-improvement result.
                        {
                            const auto insert_start=Clock::now();
                            {
                                const auto status=scores.put(shard,std::move(d.domains),value,t);
                                if(status==SharedScores::Insert::added) ++w.entries;
                                else if(status==SharedScores::Insert::duplicate) ++w.cache_duplicates;
                                else ++w.cache_refused;
                            }
                            w.key_ms+=ms(insert_start);
                        }
                    }
                }
                out.values[id]=value;
                if(!value)++w.zero;
                const auto edge_ms=ms(edge_start);w.edge_ms+=edge_ms;++w.completed;
                
                
            }
            w.local_core_calls=domains.counters().core_calls;w.local_core_removed=domains.counters().core_removed;
            w.local_core_scans=domains.counters().adjacency_scans;
            w.sparse_requests=domains.counters().sparse_requests;
        }
    };
    std::vector<std::thread> pool;
    { for(unsigned t=0;t+1<threads;++t) pool.emplace_back(run,t); run(threads-1); }
    for(auto& x:pool)x.join();
    out.timing["dynamic_compute_wall_ms"]=ms(compute_started);
    out.counters["pattern_threads_started"]=threads;
    out.counters["pivot_memory_bytes_per_region"]=pivot_memory_bytes;
    for(unsigned t=0;t<threads;++t) {
        const auto& w=workers[t];
        out.counters["worker_"+std::to_string(t)+"_tasks"]=w.completed.load();
        out.timing["worker_"+std::to_string(t)+"_edge_sum_ms"]=w.edge_ms;
        accumulate_backend_counters(out.counters,w.solver_stats);
        out.counters["pattern_solver_calls"]+=w.calls;out.counters["pattern_solver_nodes"]+=w.nodes;
        out.counters["pattern_cache_hits"]+=w.hits;out.counters["pattern_cache_entries"]+=w.entries;
        out.counters["pattern_cache_refusals"]+=w.cache_refused;
        {
            out.counters["dynamic_shared_cache_lookups"]+=w.cache_lookups;
            out.counters["dynamic_shared_cache_hits"]+=w.hits;
            out.counters["dynamic_shared_cache_misses"]+=w.cache_lookups-w.hits;
            out.counters["dynamic_shared_cache_inserts"]+=w.entries;
            out.counters["dynamic_shared_cache_duplicate_inserts"]+=w.cache_duplicates;
            out.counters["dynamic_shared_cache_refusals"]+=w.cache_refused;
            out.counters["dynamic_shared_cache_cross_worker_hits"]+=w.cache_foreign_hits;
        }
        out.counters["pattern_certified_edges"]+=w.certified;
        out.counters["pattern_lower_positive_edges"]+=w.positive;
        out.counters["pattern_warm_cutoff_calls"]+=w.warm;
        out.counters["pattern_warm_cutoff_positive_calls"]+=w.warm_positive;
        out.counters["pattern_warm_cutoff_no_improvement"]+=w.no_improvement;
        out.counters["pattern_published_witnesses"]+=w.published;
        out.counters["pattern_publications"]+=w.publications;
        out.counters["pattern_improvement_core_calls"]+=w.improvement_core;
        out.counters["pattern_improvement_raised_edges"]+=w.improvement_raised;
        out.counters["pattern_improvement_saved_solver_calls"]+=w.improvement_rejected;
        out.counters["dynamic_score_cache_bytes"]+=w.cache_bytes;
        out.counters["dynamic_local_core_calls"]+=w.local_core_calls;
        out.counters["dynamic_local_core_removed"]+=w.local_core_removed;
        out.counters["dynamic_local_core_scans"]+=w.local_core_scans;
        out.counters["domain_sparse_requests"]+=w.sparse_requests;
        out.counters["dynamic_local_partner_work"]+=w.partner.work;
        out.counters["dynamic_local_partner_removed"]+=w.partner.removed;
        out.counters["dynamic_local_partner_core_removed"]+=w.partner.core_removed;
        out.counters["dynamic_local_partner_budget_stops"]+=w.partner.budget_stops;
        out.counters["pattern_zero_edges"]+=w.zero;
        out.timing["dynamic_edge_sum_ms"]+=w.edge_ms;out.timing["dynamic_solver_sum_ms"]+=w.solver_ms;
        out.timing["dynamic_partner_sum_ms"]+=w.partner_ms;out.timing["dynamic_domains_sum_ms"]+=w.domain_ms;
        out.timing["dynamic_publish_sum_ms"]+=w.publish_ms;
        out.timing["dynamic_improvement_sum_ms"]+=w.improvement_ms;
        out.timing["dynamic_cache_sum_ms"]+=w.key_ms;
    }
    out.counters["dynamic_shared_cache_entries"]=out.counters["pattern_cache_entries"];
    out.counters["dynamic_shared_cache_bytes"]=scores.bytes();
    out.counters["dynamic_score_cache_bytes"]=scores.bytes();
    out.counters["pattern_zero_edges"]+=out.counters["pattern_global_core_skipped"];
    out.logical_payload_bytes=out.counters["dynamic_shared_bytes"]+out.counters["dynamic_score_cache_bytes"]+
        out.values.capacity()*sizeof(Score)+ids.capacity()*sizeof(std::size_t)+lower.capacity()*sizeof(std::atomic<Score>);
    out.timing["dynamic_total_ms"]=ms(started);
    out.timing["dynamic_cpu_ms"]=1000.0*(std::clock()-cpu_started)/CLOCKS_PER_SEC;
    return out;
}
} // namespace mcviz
