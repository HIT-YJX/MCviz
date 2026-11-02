#include "query_seed_heuristic.h"
#include <algorithm>
#include <limits>

namespace mcviz {
namespace {
using Clock = std::chrono::steady_clock;
struct Candidate { std::size_t node, size; Score height; };

// Build feasible seeds by sparse neighbor intersections and quota-aware choices.
// Only complete, validated witnesses raise the later exact-search cutoff.
// Work sums sampled sparse intersections; storage is local domains and bounded DFS frames.
class Greedy {
    const Graph& g;
    const PatternIndex& index;
    const QueryRequest& q;
    const Required& fixed;
    const Sets& core;
    const LayerPattern& p;
    Clock::time_point deadline;
    SearchResult& best;
    RunInfo& run;
    std::uint64_t work = 0;
    bool stopped = false;

    // Check example: an expired deadline stops at the first forced checkpoint;
    // no partial candidate is published and the exact search keeps its coverage.
    bool check(bool force = false) {
        if ((++work % 256 == 0 || force) && Clock::now() >= deadline) stopped = true;
        return !stopped;
    }
    bool has(const std::vector<VertexPos>& d, VertexPos u) const {
        return std::binary_search(d.begin(), d.end(), u);
    }
    std::vector<VertexPos> common(const std::vector<VertexPos>& a, const std::vector<VertexPos>& b) {
        const auto& small = a.size() < b.size() ? a : b;
        const auto& large = a.size() < b.size() ? b : a;
        std::vector<VertexPos> c;
        for (auto u : small) {
            if (!check()) return {};
            if (has(large,u)) c.push_back(u);
        }
        return c;
    }
    Sets copy(const Sets& d) {
        Sets out(g.K);
        for (Layer i=0;i<g.K;++i) for (auto u:d[i]) {
            if (!check()) return {};
            out[i].push_back(u);
        }
        return out;
    }
    bool select(Sets& d, Sets& s, Layer i, VertexPos u) {
        if (!has(d[i],u)) return false;
        if (!has(s[i],u)) s[i].insert(std::lower_bound(s[i].begin(),s[i].end(),u),u);
        for (auto j:p.neighbors[i]) {
            d[j]=common(d[j],g.neighbors(i,u,j));
            if (stopped) return false;
            if (d[j].size()<q.gamma[j]) return false;
            for (auto v:s[j]) {
                if (!check()) return false;
                if(!has(d[j],v)) return false;
            }
        }
        return true;
    }
    bool offer(Sets s) {
        if (!check(true)) return false;
        if (!enough_capacity(s,q.gamma)||!contains_required(s,fixed)) return false;
        for (auto [i,j]:p.edges) for (auto u:s[i]) for (auto v:s[j]) {
            if (!check()) return false;
            if(!g.has_edge(i,u,j,v)) return false;
        }
        auto value=capacity_bound(p,index.identity.spec.objective,s);
        if (!check(true)) return false;
        if (!best.found||value>best.score) {
            best={true,value,std::move(s),Status::incomplete,{}};
            ++run.counters["query_heuristic_improvements"];
            return true;
        }
        return false;
    }
    void trial(const Sets& domain, unsigned turn) {
        if (!check(true)) return;
        if (!enough_capacity(domain,q.gamma)||!contains_required(domain,fixed)) return;
        ++run.counters["query_heuristic_trials"];
        auto d=copy(domain);
        if (stopped) return;
        Sets s(g.K);
        for(auto v:fixed) if(!select(d,s,v.layer,v.pos)) return;
        unsigned step=0;
        for (;;) {
            if (!check(true)) return;
            Layer i=g.K;
            double scarcity=std::numeric_limits<double>::infinity();
            for(Layer j=0;j<g.K;++j) if(s[j].size()<q.gamma[j]) {
                double ratio=double(d[j].size()-s[j].size())/double(q.gamma[j]-s[j].size());
                if(ratio<scarcity) {scarcity=ratio;i=j;}
            }
            if(i==g.K) break;
            struct Choice {VertexPos u;double minimum,total;};
            std::vector<Choice> choices;
            const auto count=std::min<std::size_t>(128,d[i].size());
            for(std::size_t t=0;t<count;++t) {
                if (!check()) return;
                auto pos=d[i].size()<=128 ? t : (t*(d[i].size()/count)+turn*17)%d[i].size();
                auto u=d[i][pos];
                if(has(s[i],u)) continue;
                double minimum=std::numeric_limits<double>::infinity(),total=0;
                bool valid=true;
                for(auto j:p.neighbors[i]) {
                    auto c=common(d[j],g.neighbors(i,u,j));
                    if (stopped) return;
                    if(c.size()<q.gamma[j]) {valid=false;break;}
                    for(auto v:s[j]) if(!has(c,v)) {valid=false;break;}
                    double ratio=double(c.size())/double(q.gamma[j]);
                    minimum=std::min(minimum,ratio);total+=ratio;
                }
                if(valid) choices.push_back({u,minimum,total});
            }
            if(choices.empty()) return;
            std::sort(choices.begin(),choices.end(),[](const Choice& a,const Choice& b) {
                if(a.minimum!=b.minimum) return a.minimum>b.minimum;
                if(a.total!=b.total) return a.total>b.total;
                return a.u<b.u;
            });
            auto rank=step++==0 ? turn%std::min<std::size_t>(8,choices.size()) : 0;
            if(!select(d,s,i,choices[rank].u)) return;
        }
        complete(std::move(d),std::move(s),turn);
    }
    void complete(Sets d,Sets s,unsigned turn) {
        auto seed=copy(s);
        if (stopped) return;
        bool improved=offer(std::move(seed));
        if (stopped) return;
        if(improved) {
            d=copy(core);
            if (stopped) return;
            for(Layer i=0;i<g.K;++i) for(auto u:s[i])
                for(auto j:p.neighbors[i]) {
                    d[j]=common(d[j],g.neighbors(i,u,j));
                    if (stopped) return;
                }
        }
        for(Layer k=0;k<g.K;++k) {
            Layer i=(k+turn)%g.K;
            std::vector<VertexPos> candidates;
            for(auto u:d[i]) {
                if (!check()) return;
                candidates.push_back(u);
            }
            for(auto u:candidates) {
                if (!check()) return;
                if(!has(s[i],u)&&!select(d,s,i,u)) return;
            }
            auto candidate=copy(s);
            if (stopped) return;
            offer(std::move(candidate));
            if (stopped) return;
        }
    }
    Sets members(std::size_t node) {
        Sets d(g.K);
        std::vector<std::size_t> todo{node};
        const auto& f=index.forest;
        while(!todo.empty()) {
            if (!check()) return {};
            auto u=todo.back();todo.pop_back();
            if(u<f.vertices.size()) {
                auto v=f.vertices[u];
                if(has(core[v.layer],v.pos)) d[v.layer].push_back(v.pos);
            } else for(auto v:f.nodes[u].children) {
                if (!check()) return {};
                todo.push_back(v);
            }
        }
        for(auto& a:d) {
            if (!check(true)) return {};
            std::sort(a.begin(),a.end());
        }
        return d;
    }
    bool eligible(std::size_t node) {
        if (!check(true)) return false;
        ++run.counters["query_heuristic_candidates_checked"];
        std::vector<Score> counts(g.K,0);
        std::vector<bool> present(fixed.size(),false);
        std::size_t layers=0,required=0;
        std::vector<std::size_t> todo{node};
        const auto& f=index.forest;
        while(!todo.empty()) {
            if (!check()) return false;
            auto u=todo.back();todo.pop_back();
            if(u<f.vertices.size()) {
                auto v=f.vertices[u];
                if(!has(core[v.layer],v.pos)) continue;
                if(counts[v.layer]<q.gamma[v.layer]&&++counts[v.layer]==q.gamma[v.layer]) ++layers;
                for(std::size_t k=0;k<fixed.size();++k) {
                    if (!check()) return false;
                    if(!present[k]&&fixed[k].layer==v.layer&&fixed[k].pos==v.pos) {present[k]=true;++required;}
                }
                if(layers==g.K&&required==fixed.size()) return true;
            } else for(auto v:f.nodes[u].children) {
                if (!check()) return false;
                todo.push_back(v);
            }
        }
        ++run.counters["query_heuristic_candidates_rejected"];
        return false;
    }
    Layer branch_layer(const Sets& d,const Sets& s) {
        Layer i=g.K;
        double slack=std::numeric_limits<double>::infinity();
        for(Layer j=0;j<g.K;++j) if(s[j].size()<q.gamma[j]) {
            if (!check()) return g.K;
            double ratio=double(d[j].size()-s[j].size())/double(q.gamma[j]-s[j].size());
            if(i==g.K||p.neighbors[j].size()>p.neighbors[i].size()||
                (p.neighbors[j].size()==p.neighbors[i].size()&&ratio<slack)) {i=j;slack=ratio;}
        }
        return i;
    }
    std::vector<VertexPos> alternatives(const Sets& d,const Sets& s,Layer i) {
        struct Choice {VertexPos u;double minimum,total;};
        std::vector<Choice> choices;
        auto better=[](const Choice& a,const Choice& b) {
            if(a.minimum!=b.minimum) return a.minimum>b.minimum;
            if(a.total!=b.total) return a.total>b.total;
            return a.u<b.u;
        };
        const auto count=d[i].size()<=4096 ? d[i].size() : 1024;
        for(std::size_t t=0;t<count;++t) {
            if (!check()) return {};
            auto u=d[i][d[i].size()<=4096 ? t : t*(d[i].size()/count)];
            if(has(s[i],u)) continue;
            double minimum=std::numeric_limits<double>::infinity(),total=0;
            bool valid=true;
            for(auto j:p.neighbors[i]) {
                auto c=common(d[j],g.neighbors(i,u,j));
                if (stopped) return {};
                if(c.size()<q.gamma[j]) {valid=false;break;}
                for(auto v:s[j]) {
                    if (!check()) return {};
                    if(!has(c,v)) {valid=false;break;}
                }
                double ratio=double(c.size())/double(q.gamma[j]);
                minimum=std::min(minimum,ratio);total+=ratio;
            }
            if(!valid) continue;
            choices.push_back({u,minimum,total});
            std::sort(choices.begin(),choices.end(),better);
            if(choices.size()>8) choices.pop_back();
        }
        std::vector<VertexPos> out;
        for(auto c:choices) out.push_back(c.u);
        return out;
    }
    // Explore at most 4096 sparse states on a 256-choice path with 64 MiB admission.
    // Work sums intersection costs over admitted states; refusal leaves exact coverage intact.
    void backtrack(const Sets& domain) {
        if (!check(true)) return;
        if(best.found||!enough_capacity(domain,q.gamma)||!contains_required(domain,fixed)) return;
        std::uint64_t need=0;
        for(auto n:q.gamma) {
            if(n>256-need) {++run.counters["query_heuristic_backtrack_depth_skips"];return;}
            need+=n;
        }
        struct Frame {
            Sets d,s;
            std::vector<VertexPos> choices;
            std::size_t next=0,bytes=0;
            Layer layer=0;
            bool prepared=false;
        };
        constexpr std::size_t limit=64ULL<<20;
        auto estimate=[&](const Sets& d,const Sets& s) {
            std::size_t bytes=2*sizeof(Frame)+64+2*g.K*sizeof(std::vector<VertexPos>);
            for(const auto& a:d) {
                if (!check()) return std::size_t(0);
                bytes+=2*a.size()*sizeof(VertexPos);
            }
            for(const auto& a:s) {
                if (!check()) return std::size_t(0);
                bytes+=2*a.size()*sizeof(VertexPos);
            }
            return bytes;
        };
        Sets selected(g.K);
        auto root_bytes=estimate(domain,selected);
        if (stopped) return;
        std::size_t scratch=128ULL<<10;
        for(const auto& a:core) {
            if (!check()) return;
            scratch+=16*a.size()*sizeof(VertexPos);
        }
        if(scratch>=limit||root_bytes>limit-scratch) {
            ++run.counters["query_heuristic_backtrack_memory_skips"];return;
        }
        ++run.counters["query_heuristic_backtrack_calls"];
        auto root=copy(domain);
        if (stopped) return;
        for(auto v:fixed) if(!select(root,selected,v.layer,v.pos)) return;
        std::vector<Frame> stack;
        stack.reserve(257);
        stack.push_back({std::move(root),std::move(selected),{},0,root_bytes,0,false});
        std::size_t bytes=scratch+root_bytes;
        auto& states=run.counters["query_heuristic_backtrack_states"];
        while(!stack.empty()&&states<4096) {
            if (!check(true)) return;
            run.counters["query_heuristic_backtrack_peak_bytes"]=std::max<std::uint64_t>(
                run.counters["query_heuristic_backtrack_peak_bytes"],bytes);
            auto& frame=stack.back();
            if(!frame.prepared) {
                ++states;
                frame.layer=branch_layer(frame.d,frame.s);
                if (stopped) return;
                if(frame.layer==g.K) {
                    auto d=std::move(frame.d),s=std::move(frame.s);
                    stack.clear();
                    complete(std::move(d),std::move(s),0);
                    return;
                }
                frame.choices=alternatives(frame.d,frame.s,frame.layer);
                if (stopped) return;
                frame.prepared=true;
            }
            if(frame.next==frame.choices.size()) {bytes-=frame.bytes;stack.pop_back();continue;}
            auto u=frame.choices[frame.next++];
            auto child_bytes=estimate(frame.d,frame.s)+2*sizeof(VertexPos);
            if (stopped) return;
            if(child_bytes>limit-bytes) {
                ++run.counters["query_heuristic_backtrack_memory_skips"];return;
            }
            auto d=copy(frame.d);
            if (stopped) return;
            auto s=copy(frame.s);
            if (stopped) return;
            if(!select(d,s,frame.layer,u)) {
                if (stopped) return;
                continue;
            }
            ++run.counters["query_heuristic_backtrack_branches"];
            bytes+=child_bytes;
            stack.push_back({std::move(d),std::move(s),{},0,child_bytes,0,false});
        }
        if(states>=4096) run.counters["query_heuristic_backtrack_state_limit"]=1;
    }
public:
    Greedy(const Graph& graph,const PatternIndex& idx,const QueryRequest& query,const Required& required,
        const Sets& domains,Clock::time_point until,SearchResult& result,RunInfo& info):
        g(graph),index(idx),q(query),fixed(required),core(domains),p(idx.identity.spec.pattern),
        deadline(until),best(result),run(info) {}

    // Rank at most 64 eligible subtrees, then try their domains and the initial core.
    // Preparation visits at most 65536 hierarchy links; all phases share one deadline.
    void run_trials() {
        if (!check(true)) return;
        const auto& f=index.forest;
        const auto leaves=f.vertices.size();
        std::vector<std::size_t> sizes;
        std::vector<Candidate> candidates;
        std::uint64_t needed=0;
        for(auto n:q.gamma) {if(n>4096-needed) {needed=4097;break;}needed+=n;}
        auto better=[](const Candidate& a,const Candidate& b) {
            if(a.height!=b.height) return a.height>b.height;
            return a.size!=b.size ? a.size<b.size : a.node<b.node;
        };
        std::uint64_t scanned=0;
        for(std::size_t u=leaves;u<f.nodes.size()&&scanned<65536;++u) {
            if (!check()) return;
            std::size_t size=0;
            for(auto v:f.nodes[u].children) {
                if (!check()) return;
                ++scanned;size=std::min<std::size_t>(4097,size+(v<leaves?1:sizes[v-leaves]));
            }
            sizes.push_back(size);
            if(size>=needed&&size<=4096) {
                Candidate candidate{u,size,f.nodes[u].height};
                if(candidates.size()==64&&!better(candidate,candidates.back())) continue;
                if(!eligible(u)) {
                    if (stopped) return;
                    continue;
                }
                candidates.push_back(candidate);
                std::sort(candidates.begin(),candidates.end(),better);
                if(candidates.size()>64) candidates.pop_back();
            }
        }
        run.counters["query_heuristic_scanned"]=scanned;
        run.counters["query_heuristic_candidates"]=candidates.size();
        for(const auto& c:candidates) {
            if (!check(true)) return;
            auto d=members(c.node);
            if (stopped) return;
            for(unsigned t=0;t<8;++t) {
                trial(d,t);
                if (stopped) return;
            }
        }
        for(unsigned t=0;t<64;++t) {
            trial(core,t);
            if (stopped) return;
        }
        if(!best.found) {
            backtrack(core);
            if (stopped) return;
            for(const auto& c:candidates) {
                if(best.found||run.counters["query_heuristic_backtrack_states"]>=4096) break;
                if (!check(true)) return;
                auto d=members(c.node);
                if (stopped) return;
                backtrack(d);
                if (stopped) return;
            }
        }
    }
    void execute() {
        run_trials();
        if (stopped) run.counters["query_heuristic_budget_exhausted"]=1;
        run.counters["query_heuristic_work"]=work;
    }
};
}

void seed_core_heuristic(const Graph& g,const PatternIndex& index,const QueryRequest& q,const Required& fixed,
    const Sets& core,Clock::time_point deadline,SearchResult& best,RunInfo& run) {
    Greedy search(g,index,q,fixed,core,deadline,best,run);
    search.execute();
}
} // namespace mcviz
