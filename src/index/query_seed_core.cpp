#include "query_seed_core.h"
#include "query_seed_heuristic.h"
#include <algorithm>

namespace mcviz {
namespace {
using Clock = std::chrono::steady_clock;
double elapsed(Clock::time_point t) {return std::chrono::duration<double,std::milli>(Clock::now()-t).count();}

// Every improving connected witness lies in one threshold component by PA-11.
// Local gamma peeling preserves that witness by the first-deleted-vertex argument.
// Frozen domains are disjoint; the seed remains outside this strict-improvement
// search, so a vanished final graph still returns the validated incumbent.
class SeedCore {
    const Graph& g;
    const PatternIndex& index;
    const QueryRequest& q;
    const Required& fixed;
    const LayerPattern& p;
    SolverOptions options;
    SearchResult& best;
    SolverStats& stats;
    RunInfo& run;
    CoreWorkspace membership, scratch;
    std::vector<Sets> regions;
    Score minimum=0;

    bool expired() const {return options.deadline&&Clock::now()>=*options.deadline;}
    // A threshold cut materializes disjoint subtrees, intersects the prior core,
    // then peels each domain on ORIGINAL P adjacency. O(H+N) navigation plus
    // queue-core adjacency scans, O(N+H+sum |D|) live memory; no saved weight edits.
    bool collect(Score threshold,bool strict,const std::string& stage) {
        regions.clear();
        auto todo=index.forest.roots;
        const auto& f=index.forest;
        while(!todo.empty()) {
            if(expired()) return false;
            auto u=todo.back();todo.pop_back();
            if(u<f.vertices.size()) continue;
            auto height=f.nodes[u].height;
            if(strict ? height<=threshold : height<threshold) {
                for(auto v:f.nodes[u].children) todo.push_back(v);
                continue;
            }
            auto d=f.members(u,g.K);
            for(Layer i=0;i<g.K;++i) d[i].erase(std::remove_if(d[i].begin(),d[i].end(),[&](VertexPos v) {
                return !membership.find(i,v);
            }),d[i].end());
            if(!enough_capacity(d,q.gamma)||!contains_required(d,fixed)) continue;
            auto core=reduce(g,p,q.gamma,std::move(d),fixed,CoreBackend::queue,&scratch);
            ++run.counters["query_"+stage+"_core_calls"];
            run.counters["query_"+stage+"_core_removed"]+=core.removed;
            if(expired()) return false;
            if(!core.possible) continue;
            regions.push_back(std::move(core.domains));
        }
        return !expired();
    }
    // The initial core union is only a membership filter; final cutting still
    // uses the original hierarchy, preserving its proven edge-weight coverage.
    Sets unite() {
        Sets d(g.K);
        for(auto& region:regions) for(Layer i=0;i<g.K;++i)
            d[i].insert(d[i].end(),region[i].begin(),region[i].end());
        for(auto& a:d) std::sort(a.begin(),a.end());
        return d;
    }
public:
    SeedCore(const Graph& graph,const PatternIndex& idx,const QueryRequest& query,const Required& required,
        SolverOptions search,SearchResult& result,SolverStats& statistics,RunInfo& info):
        g(graph),index(idx),q(query),fixed(required),p(idx.identity.spec.pattern),options(search),
        best(result),stats(statistics),run(info) {}

    /* Regression example: a zero heuristic budget skips seed construction.
       Both cuts still cover every improving clique, so the exact answer is unchanged. */
    bool execute() {
        run.counters["query_seed_core_used"]=1;
        run.counters["query_seed_core_queue_core"]=1;
        run.counters["query_heuristic_budget_ms"]=index.options.query_heuristic_ms;
        for(auto key:{"query_seed_core_initial_found","query_seed_core_initial_score","query_heuristic_trials",
                      "query_heuristic_improvements","query_heuristic_budget_exhausted",
                      "query_heuristic_candidates_checked","query_heuristic_candidates_rejected"}) run.counters[key]=0;
        for(auto stage:{"initial_core","final_graph"}) {
            run.counters[std::string("query_")+stage+"_complete"]=0;
            run.counters[std::string("query_")+stage+"_sampled"]=0;
        }
        for(auto key:{"query_initial_hierarchy_core_ms","query_heuristic_ms","query_final_cut_core_ms",
                      "query_preprocessing_total_ms","query_exact_search_ms"}) run.timing[key]=0;
        auto start=Clock::now();
        membership.bind(g,index.core_domains);
        bool possible=enough_capacity(index.core_domains,q.gamma)&&contains_required(index.core_domains,fixed);
        bool complete=true;
        if(possible) {
            minimum=score(p,index.identity.spec.objective,q.gamma);
            complete=collect(minimum,false,"initial_core");
        }
        Sets core=unite();
        membership.bind(g,core);
        run.counters["query_initial_core_complete"]=complete;
        run.counters["query_initial_core_threshold"]=minimum;
        run.timing["query_initial_hierarchy_core_ms"]=elapsed(start);
        if(!complete||expired()) return false;

        start=Clock::now();
        if(index.options.query_heuristic_ms&&!regions.empty()) {
            auto deadline=start+std::chrono::milliseconds(index.options.query_heuristic_ms);
            if(options.deadline) deadline=std::min(deadline,*options.deadline);
            seed_core_heuristic(g,index,q,fixed,core,deadline,best,run);
        }
        run.timing["query_heuristic_ms"]=elapsed(start);
        run.counters["query_seed_core_initial_found"]=best.found;
        run.counters["query_seed_core_initial_score"]=best.found?best.score:0;
        run.counters["query_main_initial_score"]=best.found?best.score:0;
        run.counters["query_final_graph_threshold"]=best.found?best.score:minimum;
        run.counters["query_final_graph_strict"]=best.found;
        if(expired()) return false;

        start=Clock::now();
        complete=possible ? collect(best.found?best.score:minimum,best.found,"final_graph") : true;
        run.counters["query_final_graph_complete"]=complete;
        run.timing["query_final_cut_core_ms"]=elapsed(start);
        if(!complete||expired()) return false;

        // The full final graph and its accounting are frozen at this point.
        // Later incumbents can prune DFS but never rewrite the measured graph.
        start=Clock::now();
        for(auto& d:regions) {
            if(expired()||(options.node_limit&&stats.nodes>=*options.node_limit)) {complete=false;break;}
            auto search=options;
            search.core_backend=CoreBackend::queue;
            search.core_workspace=&scratch;
            search.research_skip_initial_core=true;
            if(search.node_limit) search.node_limit=*search.node_limit-stats.nodes;
            auto r=solve_region(g,p,index.identity.spec.objective,q.gamma,std::move(d),fixed,
                                best.found?std::optional<Score>(best.score):std::nullopt,search);
            ++run.solver_calls;++run.regions;++run.counters["query_main_solver_calls"];
            run.solver=run.solver=="none"||run.solver==r.backend ? r.backend : "mixed";
            merge_solver_stats(stats,r.stats);
            if(r.best.found&&(!best.found||r.best.score>best.score)) best=std::move(r.best);
            if(!r.complete) {complete=false;break;}
        }
        run.timing["query_exact_search_ms"]=elapsed(start);
        return complete;
    }
};
}

bool seed_core_query(const Graph& g,const PatternIndex& index,const QueryRequest& q,const Required& fixed,
    SolverOptions options,SearchResult& best,SolverStats& stats,RunInfo& run) {
    SeedCore search(g,index,q,fixed,options,best,stats,run);
    bool complete=search.execute();
    run.timing["query_preprocessing_total_ms"]=run.timing["query_initial_hierarchy_core_ms"]+
        run.timing["query_heuristic_ms"]+run.timing["query_final_cut_core_ms"];
    run.timing["query_main_ms"]=run.timing["query_exact_search_ms"];
    run.timing["query_seed_ms"]=run.timing["query_heuristic_ms"];
    run.counters["query_main_nodes"]=stats.nodes;
    run.counters["query_main_pivot_nodes"]=stats.pivot_nodes;
    run.counters["query_main_cover_nodes"]=stats.nodes-stats.pivot_nodes;
    run.counters["query_total_nodes"]=stats.nodes;
    return complete;
}
} // namespace mcviz
