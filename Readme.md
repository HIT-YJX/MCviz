# MCviz

MCviz finds maximum layer-pattern cliques in static multipartite graphs. This C++17 release provides MCviz, MCviz+, BCviz-M, Motif-LP, KPC-LP, and KSPC-LP, with index construction, hierarchy conversion, and exact search. WordNet and LastFM are bundled; [datasets/Readme.md](datasets/Readme.md) describes the eight paper datasets and their sources.

## Problem and algorithms

A fixed layer pattern `P` specifies the layer pairs that require complete connections. An answer selects at least `gamma[i]` vertices in every layer and connects every selected pair across each edge of `P`. Other layer pairs impose no constraint. The `edge` objective maximizes the number of selected pattern edges; `vertex` maximizes the total selected vertices; and `balanced` maximizes the size of the smallest selected layer. All release algorithms support these three objectives. Thresholds are positive, and indexed queries require `gamma[i] >= beta[i]`, where `beta` is the construction threshold vector.

MCviz assigns each pattern edge its exact cohesion, the maximum feasible objective score of a clique containing the edge. MCviz+ uses safe upper bounds in place of exact cohesion. Both build a maximum spanning forest and a component hierarchy to locate candidate regions. Their default query first cuts the hierarchy and applies local threshold cores, finds an initial feasible clique with a bounded heuristic, then cuts again to retain possible strict improvements. Exact regional search completes the answer; exhausting the heuristic budget does not discard unresolved search coverage.

| Method | Construction or search principle |
| --- | --- |
| MCviz | Exact cohesion: dynamic construction for `edge`; shared and anchor-based star enumeration for `vertex` and `balanced`, with exact per-edge fallback. |
| MCviz+ | Automatic upper bounds: support bounds for stars, structural bounds for trees, and safe bounds for other supported patterns. |
| BCviz-M | Degree-based edge bounds and a global vertex order organize candidate regions for exact search. |
| Motif-LP | Candidate/exclusion enumeration with pivots retains the maximum feasible compatible set. |
| KPC-LP | Ordered per-layer seeds and extensions use threshold and score-capacity pruning. |
| KSPC-LP | Star queries combine maximal bicliques from center--leaf pairs; other supported queries use exhaustive subset search. |

The four baselines are the project's layer-pattern adaptations. KSPC-LP requires all layers through `--min-layers K`. Search can take exponential time. Configured work and memory limits retain their algorithmic fallbacks, including uncached computation and cover search when Pivot admission fails. The release omits experimental strategy switches and exception recovery.

Patterns are `star`, `complete`, `path`, `cycle`, `empty`, or `file:path` with one 0-based layer pair per line. `--primary` defines the star center; for other patterns it is only an optional anchor. `--required layer:id` requests a clique containing that vertex. The connected-pattern index guarantees apply to connected patterns; disconnected and empty patterns use separate fallback paths.

## Structure and build

All algorithm implementations are C++. CMake 3.16 or newer and a C++17 compiler with standard thread support build the executables below.

```text
src/common/       Graph input, patterns, objectives, CLI, and timing tables
src/search/       Exact regional search and construction helpers
src/index/        MCviz/MCviz+ indexes, queries, hierarchy, and conversion
src/mcviz/        MCviz entry point
src/mcviz_plus/   MCviz+ entry point
src/baseline/     BCviz-M, Motif-LP, KPC-LP, and KSPC-LP
datasets/         WordNet and LastFM inputs
Index/            Generated and reused indexes
log/              Timing TXT tables
```

From this directory, run:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 1
```

Executables are written to `build/bin/`. A graph or vertex-table basename is resolved under `datasets/` when it is absent from the working directory. The default index and log directories are `Index/` and `log/` under the release source directory used by CMake. The commands below use these default directories directly, without run-specific subdirectories. Run each command separately; timing records append to the corresponding TXT file in `log/`. No separate cohesion files are written.

## Defaults and execution

The exported policies follow the production defaults of November 1, 2026. MCviz exact `edge` construction uses `general-dynamic-all` with 32 workers on connected nonempty patterns. For `vertex` and `balanced`, star construction uses shared and anchor-based enumeration with exact per-edge fallback when configured budgets are exhausted; other patterns retain exact regional construction. These non-edge construction routes use one worker. MCviz+ requests 32 workers for applicable star-support construction; nonparallel upper-bound routes use one worker. MCviz and MCviz+ queries use `seed-core`, automatic Pivot with dynamic ordering, a 512 MiB per-region Pivot budget, and a 1,000 ms initial-solution budget. Search is single-threaded. Actual worker counts may be limited by task size and cache admission.

Use `--build-threads 1` for small local runs. `--pivot-memory-mb` and `--query-heuristic-ms` set the search budgets; a zero Pivot budget selects cover fallback. The paper's larger Freebase, DRKG, and AMiner heuristic budgets require explicit `--query-heuristic-ms 60000`; the program does not choose them from dataset names.

The following Bash example builds and queries WordNet-P3. Its explicit vertex table preserves the complete projected domain, including isolated vertices. The paper's main query threshold is `4,4,4`; use `2,2,2` to run the smaller-threshold example on the same indexes.

```bash
args=(--graph WordNet-P3.txt --vertices WordNet-P3.vertices.tsv
      --layers 3 --pattern path --primary 1 --objective edge)

build/bin/mcviz_layer construct "${args[@]}" --beta 2,2,2 --build-threads 1
build/bin/mcviz_layer query "${args[@]}" --beta 2,2,2 --gamma 4,4,4

build/bin/mcviz_plus_layer construct "${args[@]}" --beta 2,2,2 --build-threads 1
build/bin/mcviz_plus_layer query "${args[@]}" --beta 2,2,2 --gamma 4,4,4

# Query the same WordNet index with a smaller feasible threshold.
build/bin/mcviz_layer query "${args[@]}" --beta 2,2,2 --gamma 2,2,2

# Exercise budget fallback without the heuristic or Pivot workspace.
build/bin/mcviz_layer query "${args[@]}" --beta 2,2,2 --gamma 2,2,2 \
    --query-heuristic-ms 0 --pivot-memory-mb 0

build/bin/baseline_layer construct "${args[@]}" --method bcvizm --beta 2,2,2
build/bin/baseline_layer query "${args[@]}" --method bcvizm --beta 2,2,2 --gamma 4,4,4

for method in motif kpc; do
    build/bin/baseline_layer solve "${args[@]}" --method "$method" --gamma 4,4,4
done
build/bin/baseline_layer solve "${args[@]}" --method kspc --min-layers 3 --gamma 4,4,4
```

For MCviz+ on LastFM, run the following commands individually to build at beta=2 and query at gamma=4:

```bash
graph=(--graph LastFM.txt --layers 3 --pattern star --primary 0)

build/bin/mcviz_plus_layer construct "${graph[@]}" --objective edge --beta 2,2,2 --build-threads 1
build/bin/mcviz_plus_layer query "${graph[@]}" --objective edge --beta 2,2,2 --gamma 4,4,4

build/bin/mcviz_plus_layer construct "${graph[@]}" --objective vertex --beta 2,2,2 --build-threads 1
build/bin/mcviz_plus_layer query "${graph[@]}" --objective vertex --beta 2,2,2 --gamma 4,4,4

build/bin/mcviz_plus_layer construct "${graph[@]}" --objective balanced --beta 2,2,2 --build-threads 1
build/bin/mcviz_plus_layer query "${graph[@]}" --objective balanced --beta 2,2,2 --gamma 4,4,4
```

Use `mcviz_layer` for the corresponding MCviz commands. For WordNet-P3, replace `graph` with `(--graph WordNet-P3.txt --vertices WordNet-P3.vertices.tsv --layers 3 --pattern path --primary 1)`.

Each objective has its own index; an index built for one objective must not be used for another.

Index lookup matches the graph, pattern, objective, and construction thresholds. An explicit `--index name.idx` selects a file in `Index/`. BCviz-M also supports `solve`, which constructs its order within the call. For LastFM, replace the shared arguments and repeat the commands:

```bash
args=(--graph LastFM.txt --layers 3
      --pattern star --primary 0 --objective edge)
```

## Hierarchy conversion and output

The conversion program loads an index, converts its component hierarchy to the project's shared linear representation, and restores the hierarchy in memory. The shared representation preserves cross-layer component identities; it is distinct from independently stored per-layer orders.

After building the WordNet-P3 indexes, run:

```bash
build/bin/mcviz_representation_conversion exact query \
    --graph WordNet-P3.txt --vertices WordNet-P3.vertices.tsv \
    --layers 3 --pattern path --primary 1 --objective edge --beta 2,2,2

build/bin/mcviz_representation_conversion upper query \
    --graph WordNet-P3.txt --vertices WordNet-P3.vertices.tsv \
    --layers 3 --pattern path --primary 1 --objective edge --beta 2,2,2
```

Commands append compact typed TSV tables to `log/builds.txt`, `log/searches.txt`, and `log/conversions.txt`. Each table begins with the `MCviz-test2 TSV v1` marker, column names, and column types; `\N` represents a missing value. Records include method, input configuration, status, and timing, with result scores for searches and sizes for indexes and representations. A compatible existing index is reused rather than rebuilt; this records `status=reused` without a new construction time. Infeasibility is distinguished from interruption, and a feasible zero-score answer retains `found=true`.

Build time excludes graph loading and index writing. Search time includes graph/index loading and search from CLI entry, excluding process creation and log writing. Conversion times cover the two in-memory transformations. These timing scopes should be matched explicitly when comparing to paper plots. Repeated runs append rows; select the intended configuration and repetitions before plotting.
