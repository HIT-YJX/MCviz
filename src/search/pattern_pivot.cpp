#include "pattern_pivot.h"
#include "common/text.h"
#include <algorithm>
#include <limits>
#include <memory>
#include <new>
#include <stdexcept>

namespace mcviz {
namespace {
using Word=std::uint64_t;
struct Frame { bool entered=false; std::size_t from=0; };
struct Workspace {
    std::size_t n,w;
    std::vector<LocalVertex> vertices;
    std::vector<std::size_t> offset,path;
    std::vector<Word> adj,masks;
    std::vector<Score> counts,cap;
    std::vector<Frame> frames;
    Workspace(const Sets& d,std::size_t size):n(size),w((n+63)/64),offset(d.size()+1),
        adj(n*w),masks(3*(n+1)*w),counts(d.size()),cap(d.size()),frames(n+1) {
        vertices.reserve(n);path.reserve(n);
        for(Layer i=0;i<d.size();++i) {
            offset[i]=vertices.size();

            for(auto u:d[i])vertices.push_back({i,u});
        }
        offset.back()=vertices.size();
    }


    Word* row(std::size_t u) {return adj.data()+u*w;}



    Word* mask(std::size_t depth,unsigned which) {return masks.data()+(3*depth+which)*w;}
};
bool expired(const SolverOptions& o) {return o.deadline && std::chrono::steady_clock::now()>=*o.deadline;}

void set(Word* a,std::size_t u) {a[u/64]|=Word(1)<<(u%64);}



void clear(Word* a,std::size_t u) {a[u/64]&=~(Word(1)<<(u%64));}


bool empty(const Word* a,std::size_t w) {for(std::size_t i=0;i<w;++i)if(a[i])return false;return true;}


template<class F> void each(const Word* a,std::size_t w,F f) {
    for(std::size_t i=0;i<w;++i) {auto bits=a[i];while(bits) {auto b=__builtin_ctzll(bits);f(64*i+b);bits&=bits-1;}}
}



std::uint64_t estimate(const Sets& d,std::size_t n) {
    const long double w=(static_cast<long double>(n)+63)/64;
    const long double bytes=4096+8*(4*static_cast<long double>(n)+3)*w+
        (n+1.0L)*sizeof(Frame)+static_cast<long double>(n)*(sizeof(LocalVertex)+sizeof(std::size_t)+sizeof(VertexPos))+
        (d.size()+1.0L)*(sizeof(std::size_t)+2*sizeof(Score)+sizeof(std::vector<VertexPos>));
    return bytes>=UINT64_MAX?UINT64_MAX:static_cast<std::uint64_t>(bytes);
}



void finish(RegionResult& out,std::optional<Score> cutoff) {
    if(!out.complete)out.best.status=Status::incomplete;
    else if(out.best.found)out.best.status=Status::optimal;
    else out.best.status=cutoff?Status::incomplete:Status::infeasible;
    if(out.complete && cutoff && !out.best.found)out.best.message="Region complete: no strict improvement above cutoff";
}
}

void merge_solver_stats(SolverStats& a,const SolverStats& b) {
    a.nodes+=b.nodes;a.leaves+=b.leaves;a.core_removed+=b.core_removed;a.bound_prunes+=b.bound_prunes;

    a.pa01_calls+=b.pa01_calls;a.pa04_calls+=b.pa04_calls;a.pa09_calls+=b.pa09_calls;a.pa09_prunes+=b.pa09_prunes;
    a.pa08_tree_calls+=b.pa08_tree_calls;a.pa08_core_calls+=b.pa08_core_calls;
    a.pa08_prunes+=b.pa08_prunes;a.pa08_budget_fallbacks+=b.pa08_budget_fallbacks;
    a.cover_calls+=b.cover_calls;a.pivot_calls+=b.pivot_calls;a.pivot_nodes+=b.pivot_nodes;

    a.pivot_pivots+=b.pivot_pivots;a.pivot_branches+=b.pivot_branches;

    a.pivot_prepare_bytes=std::max(a.pivot_prepare_bytes,b.pivot_prepare_bytes);
    a.pivot_memory_fallbacks+=b.pivot_memory_fallbacks;
}


void accumulate_backend_counters(std::map<std::string,std::uint64_t>& c,const SolverStats& s) {
    c["cover_calls"]+=s.cover_calls;c["pivot_calls"]+=s.pivot_calls;c["pivot_nodes"]+=s.pivot_nodes;

    c["pivot_pivots"]+=s.pivot_pivots;c["pivot_branches"]+=s.pivot_branches;
    c["pivot_memory_fallbacks"]+=s.pivot_memory_fallbacks;
    c["pivot_prepare_bytes"]=std::max(c["pivot_prepare_bytes"],s.pivot_prepare_bytes);
}



// Search the implicit compatibility graph with a pivot and dynamic branch order.
// Time is exponential in the region size; admitted matrix and DFS masks use O(n^2/64) words.
// A refused matrix budget falls back to exact cover enumeration.
RegionResult solve_pivot_region(const Graph& g,const LayerPattern& p,Objective obj,
    const Thresholds& gamma,Sets domains,const Required& fixed,std::optional<Score> cutoff,const SolverOptions& opt) {
    RegionResult out;out.backend="pivot_ordered";out.cutoff_applied=cutoff.has_value();

    out.stats.pivot_calls=1;out.best.sets.resize(g.K);

    if(expired(opt) || (opt.node_limit && !*opt.node_limit)) {
        out.complete=false;finish(out,cutoff);return out;
    }

    CoreResult core;
    if(opt.research_skip_initial_core) {
        auto d=condition_domains(g,p,std::move(domains),fixed);

        core.domains=std::move(d.domains);core.possible=d.possible&&enough_capacity(core.domains,gamma);
    } else {
        core=reduce(g,p,gamma,std::move(domains),fixed,opt.core_backend,opt.core_workspace);
        out.stats.pa01_calls=1;out.stats.core_removed=core.removed;
    }
    if(expired(opt)) {out.complete=false;finish(out,cutoff);return out;}
    if(!core.possible) {finish(out,cutoff);return out;}
    std::size_t n=0;for(const auto& d:core.domains) {
        require(d.size()<=std::numeric_limits<std::size_t>::max()-n,"Pivot domain size overflow");n+=d.size();
    }
    out.stats.pivot_prepare_bytes=estimate(core.domains,n);

    auto fallback=[&](const std::string& reason) {
        auto options=opt;options.backend=RegionBackend::cover;options.research_skip_initial_core=true;
        auto r=solve_cover_region(g,p,obj,gamma,std::move(core.domains),fixed,cutoff,options);
        ++out.stats.pivot_memory_fallbacks;merge_solver_stats(r.stats,out.stats);r.fallback=reason;return r;
    };
    if(out.stats.pivot_prepare_bytes>opt.pivot_memory_bytes ||
       out.stats.pivot_prepare_bytes>=std::numeric_limits<std::size_t>::max())return fallback("memory_budget");
    std::unique_ptr<Workspace> state;
    {
        state=std::make_unique<Workspace>(core.domains,n);
        for(Layer i=0;i<g.K;++i)out.best.sets[i].reserve(core.domains[i].size());
    }
    auto& s=*state;
    for(std::size_t u=0;u<n;++u) {
        if(expired(opt)) {out.complete=false;break;}
        const auto a=s.vertices[u];
        for(std::size_t v=u+1;v<n;++v) {
            if((v&1023)==0 && expired(opt)) {out.complete=false;break;}
            const auto b=s.vertices[v];

            if(a.layer==b.layer || !p.adjacent(a.layer,b.layer) || g.has_edge(a.layer,a.pos,b.layer,b.pos)) {
                set(s.row(u),v);set(s.row(v),u);
            }
        }
        if(!out.complete)break;
    }

    if(!out.complete) {finish(out,cutoff);return out;}

    auto C=s.mask(0,0);

    for(std::size_t u=0;u<n;++u)set(C,u);
    for(auto q:fixed) {
        const auto& d=core.domains[q.layer];auto at=std::lower_bound(d.begin(),d.end(),q.pos);

        const auto u=s.offset[q.layer]+std::size_t(at-d.begin());
        if(std::find(s.path.begin(),s.path.end(),u)!=s.path.end())continue;

        s.path.push_back(u);++s.counts[q.layer];

        for(std::size_t z=0;z<s.w;++z)C[z]&=s.row(u)[z];
    }
    std::size_t depth=0;bool done=false;

    auto pop=[&] {
        if(depth==0) {done=true;return;}
        const auto u=s.frames[depth].from;--s.counts[s.vertices[u].layer];s.path.pop_back();--depth;
        clear(s.mask(depth,0),u);set(s.mask(depth,1),u);
    };
    auto degree=[&](std::size_t u,const Word* candidates) {
        std::size_t d=0;for(std::size_t z=0;z<s.w;++z)d+=__builtin_popcountll(candidates[z]&s.row(u)[z]);return d;
    };
    while(!done) {
        if(expired(opt)) {out.complete=false;break;}

        auto& frame=s.frames[depth];C=s.mask(depth,0);auto X=s.mask(depth,1);auto branch=s.mask(depth,2);
        if(!frame.entered) {
            if(opt.node_limit && out.stats.nodes>=*opt.node_limit) {out.complete=false;break;}
            ++out.stats.nodes;++out.stats.pivot_nodes;s.cap=s.counts;
            each(C,s.w,[&](auto u){++s.cap[s.vertices[u].layer];});
            bool enough=true;for(Layer i=0;i<g.K;++i)if(s.cap[i]<gamma[i])enough=false;

            if(!enough) {++out.stats.bound_prunes;pop();continue;}
            const auto ub=score(p,obj,s.cap);
            const auto limit=out.best.found?std::optional<Score>(out.best.score):cutoff;
            if(limit && ub<=*limit) {++out.stats.bound_prunes;pop();continue;}
            if(empty(C,s.w)) {
                if(empty(X,s.w)) {
                    ++out.stats.leaves;out.best.found=true;out.best.score=ub;
                    for(auto& d:out.best.sets)d.clear();

                    for(auto u:s.path)out.best.sets[s.vertices[u].layer].push_back(s.vertices[u].pos);
                    for(auto& d:out.best.sets)std::sort(d.begin(),d.end());
                }
                pop();continue;
            }

            std::size_t pivot=n,rank=0;
            for(std::size_t z=0;z<s.w;++z) {
                if(expired(opt)) {out.complete=false;break;}
                auto bits=C[z]|X[z];while(bits) {
                    const auto u=64*z+__builtin_ctzll(bits),d=degree(u,C);bits&=bits-1;
                    if(pivot==n || d>rank) {pivot=u;rank=d;}
                }
            }

            if(!out.complete)break;

            ++out.stats.pivot_pivots;
            for(std::size_t z=0;z<s.w;++z)branch[z]=C[z]&~s.row(pivot)[z];

            frame.entered=true;
        }
        if(empty(branch,s.w)) {pop();continue;}
        std::size_t v=n,rank=0;

        for(std::size_t z=0;z<s.w;++z) {
            if(expired(opt)) {out.complete=false;break;}
            auto bits=branch[z];while(bits) {
                const auto u=64*z+__builtin_ctzll(bits),d=degree(u,C);bits&=bits-1;
                if(v==n || d>rank) {v=u;rank=d;}
            }
        }

        if(!out.complete)break;
        clear(branch,v);++out.stats.pivot_branches;
        auto nextC=s.mask(depth+1,0),nextX=s.mask(depth+1,1);

        for(std::size_t z=0;z<s.w;++z) {nextC[z]=C[z]&s.row(v)[z];nextX[z]=X[z]&s.row(v)[z];}
        s.path.push_back(v);++s.counts[s.vertices[v].layer];++depth;s.frames[depth]={false,v};
    }

    state.reset();
    finish(out,cutoff);
    return out;
}


// Dispatch regional search to the selected Pivot or cover backend.
RegionResult solve_region(const Graph& g,const LayerPattern& p,Objective obj,const Thresholds& gamma,
    Sets domains,const Required& fixed,std::optional<Score> cutoff,const SolverOptions& opt) {
    if(opt.backend==RegionBackend::pivot_ordered)
        return solve_pivot_region(g,p,obj,gamma,std::move(domains),fixed,cutoff,opt);
    return solve_cover_region(g,p,obj,gamma,std::move(domains),fixed,cutoff,opt);
}
}
