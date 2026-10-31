# MCviz

MCviz searches for maximum layer-pattern cliques in static multipartite graphs. This release contains MCviz, MCviz+, and four layer-pattern baselines, together with index construction, hierarchy conversion, and exact search. WordNet and LastFM are included; [datasets/Readme.md](datasets/Readme.md) describes all eight datasets used in the paper.

## Problem

A fixed pattern `P` specifies which layer pairs require complete connections. A feasible answer selects at least `gamma[i]` vertices from every layer; selected sets must be fully connected across each edge of `P`. Other layer pairs impose no constraint. All layer thresholds are positive.

The objective is the global maximum of one of three scores:

| Option | Score for selected sets `S[i]` |
| --- | --- |
| `edge` | Sum of `|S[i]| * |S[j]|` over pattern edges `(i,j)` |
| `vertex` | Sum of `|S[i]|` over all layers |
| `balanced` | Minimum `|S[i]|` over all layers |

Patterns may be `star`, `complete`, `path`, `cycle`, `empty`, or a file containing one 0-based layer pair per line. `--primary` chooses the center of a star; for other patterns it is an optional search anchor and does not change the constraints. `--required layer:id` restricts the answer to contain a specified vertex. Index queries require `gamma[i] >= beta[i]`, where `beta` is the construction threshold vector. The connected-pattern index and hierarchy results apply to connected `P`; disconnected and empty patterns use separate search fallback paths.

## Algorithms

MCviz and MCviz+ organize candidate regions by edge weights, then recover or search an exact answer under the original pattern constraints. A maximum spanning forest preserves threshold components, and the component hierarchy converts reversibly to a shared linear representation with linked layer intervals.

| Method | Core idea |
| --- | --- |
| MCviz | Compute each required edge's exact cohesion: the largest feasible score of a clique containing that edge. Use these weights to recover witnesses or search hierarchy regions. |
| MCviz+ | Replace exact cohesion with admissible upper bounds, including support bounds for stars. Query regions are still solved exactly. |
| BCviz-M | Build a global vertex order from safe degree-based edge bounds, then perform exact search within the resulting candidate regions. |
| Motif-LP | Enumerate maximal compatible sets with candidate/exclusion sets and pivots, then retain the best feasible score. |
| KPC-LP | Choose the smallest selected vertex of each layer as a seed, extend in vertex order, and prune using threshold and score-capacity bounds. |
| KSPC-LP | For unrestricted star queries, enumerate maximal bicliques on each center--leaf pair and combine their center intersections. Other patterns and specified-vertex queries use exhaustive subset search. |

The baselines are the layer-pattern adaptations implemented in this release. KSPC-LP requires all `K` layers (`--min-layers K`). Exact search can take exponential time. MCviz exact edge construction uses the paper's general-dynamic-all policy: queue core, partner reduction, dynamic edge scheduling, a shared completed-domain cache, and witnessed bounds. MCviz+ keeps automatic star-support or structural upper bounds. Legacy exact star/singleton builders and research strategy switches are omitted. Configured work and memory budget refusals preserve exact fallback. Actual allocation or thread creation failures terminate the process.

## Layout and build

All algorithm code is under `src/`. The build requires CMake 3.16 or later and a C++17 compiler with standard thread support.

```text
src/common/       Graph input, patterns, objectives, CLI, index files, and timing logs
src/search/       Exact regional search and construction helpers
src/index/        MCviz/MCviz+ construction, query, hierarchy, and linear conversion
src/mcviz/        MCviz executable entry
src/mcviz_plus/   MCviz+ executable entry
src/baseline/     BCviz-M, Motif-LP, KPC-LP, and KSPC-LP
datasets/         WordNet and LastFM inputs
Index/            Generated index files
log/              Generated timing tables
```

Build from the release directory:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 1
```

The executables appear in `build/bin/`. Default input, index, and log paths refer to `datasets/`, `Index/`, and `log/` beneath the source directory used when configuring CMake. A graph or vertex-table basename is resolved under `datasets/` when it does not exist in the current directory. `--index-dir` and `--log-dir` override output locations.

## Build and search

The following Bash examples use WordNet-P3 with the path `0--1--2` and the edge objective. Its explicit vertex table is required to preserve the complete projected vertex domain. Index lookup is automatic from the graph, pattern, objective, and construction thresholds; `--index name.idx` can select a file in `Index/` explicitly.

```bash
args=(--graph wordnet-p3.txt --vertices wordnet-p3.vertices.tsv
      --layers 3 --pattern path --primary 1 --objective edge)

for exe in mcviz_layer mcviz_plus_layer; do
    build/bin/$exe construct "${args[@]}" --beta 2,2,2
    build/bin/$exe query "${args[@]}" --beta 2,2,2 --gamma 2,2,2 --pivot-memory-mb 512
done

build/bin/baseline_layer construct "${args[@]}" --method bcvizm --beta 2,2,2
build/bin/baseline_layer query "${args[@]}" --method bcvizm --beta 2,2,2 --gamma 2,2,2
```

All four baselines also support an independent `solve` call. BCviz-M builds its order within that call; the other three methods search directly.

```bash
for method in bcvizm motif kpc; do
    build/bin/baseline_layer solve "${args[@]}" --method "$method" --gamma 2,2,2
done
build/bin/baseline_layer solve "${args[@]}" --method kspc --min-layers 3 --gamma 2,2,2
```

For LastFM, replace the shared arguments with the following and use the same construction/query commands:

```bash
args=(--graph lastfm_hetrec_k3.txt --layers 3
      --pattern star --primary 0 --objective edge)
```

Search uses one thread. New MCviz+ star constructions accept `--build-threads N` (default `1`, maximum `64`). MCviz and MCviz+ accept `--pivot-memory-mb N` (default `64`); a zero budget selects the existing cover-search fallback. Use `--help` on either search executable for the remaining input options.

## Hierarchy conversion and logs

The conversion executable loads a connected-pattern index, converts its hierarchy to the shared linear representation, and converts it back in memory. It records each conversion phase separately. After constructing the WordNet-P3 indexes above, run:

```bash
for mode in exact upper; do
    build/bin/mcviz_representation_conversion "$mode" query \
        --graph wordnet-p3.txt --vertices wordnet-p3.vertices.tsv \
        --layers 3 --pattern path --primary 1 --objective edge --beta 2,2,2
done
```

Commands print a short status summary and append typed TSV records to `log/builds.txt`, `log/searches.txt`, or `log/conversions.txt`. The tables use the `MCviz-test2 TSV v1` marker, followed by column names, types, and data; `\N` denotes a missing value. Threshold vectors and canonical pattern edges occupy structured cells within the TXT table.

Build time measures the construction kernel, excluding input loading and index writing. Search time starts at CLI entry and includes graph/index input and search; a BCviz-M `solve` also includes order construction. It excludes process creation and logging, so it is distinct from the paper's process-wall query timings. Conversion time measures each in-memory conversion; representation sizes are logical bytes. `value_seconds`, `value_mib`, and `value_ms` are plotting-unit aliases where provided.

Reusing an index records `status=reused` with no new construction time. A completed infeasible search records `status=ok`, `found=false`, and no score; a feasible zero-score answer retains `found=true`. Repeated runs append rows. Select the intended configuration and repetition before passing records to fixed-matrix plotters; these files contain only operations actually run, not a complete eight-dataset experiment matrix or ablation study.

## Paper construction configuration

Use a fresh index directory when measuring construction. MCviz edge builds use the following parameters from paper-exp.md; the default thread count remains one.

```bash
build/bin/mcviz_layer construct "${args[@]}" --beta 2,2,2 \
    --cohesion-strategy general-dynamic-all --solver auto --build-threads 32 \
    --build-core queue --build-filter off --cohesion-work-budget 1000000000 \
    --cohesion-cache-mb 8192 --pivot-memory-mb 512
```

For MCviz+, use automatic upper-bound construction with one worker and the default support-cache policy. Search uses `--pivot-memory-mb 512` and one thread. Historical singleton and ablation configurations in paper-exp.md are not production options in this release. Non-edge objectives retain generic exact construction; their timing is outside the paper configuration.

No exception recovery is installed. Configured cache admission refusals use sparse adjacency or uncached computation; partner-work refusal retains unresolved candidates; Pivot admission refusal uses cover search; butterfly-work refusal leaves unprocessed edges available. Actual resource allocation and thread-creation failures are fatal.
