# Datasets

This release includes WordNet and LastFM. The table uses the names and abbreviations in the paper and describes all eight processed datasets. Counts refer to the complete processed inputs in the paper's dataset table, rather than the full upstream collections or the smaller fixed-pattern projections used by some experiments.

| Dataset and source | Abbreviation | Layers | Vertices | Edges | Layer types | Paper pattern |
| --- | --- | ---: | ---: | ---: | --- | --- |
| [WordNet](https://wordnet.princeton.edu/) | WN | 25 | 29,904 | 27,544 | Synsets grouped by lexname | Path |
| [LastFM](https://grouplens.org/datasets/hetrec-2011/) | LF | 3 | 24,164 | 216,630 | User, Artist, Tag | Star |
| [Hetionet](https://het.io/) | HN | 6 | 34,134 | 1,315,516 | Gene, Anatomy, Biological process, Compound, Disease, Pathway | Star |
| [DRKG](https://github.com/gnn4dr/DRKG) | DK | 5 | 75,655 | 1,515,820 | Gene, Anatomy, Biological process, Compound, Disease | Tree |
| [BibSonomy](https://www.bibsonomy.org/) | BS | 3 | 977,914 | 3,754,828 | User, Tag, Item | Star |
| [Freebase](https://freebase-easy.cs.uni-freiburg.de/dump/) | FB | 6 | 1,669,251 | 2,536,415 | Genre, Group, Composer, Musician, Artist, Album | Star |
| [YouTube](https://snap.stanford.edu/data/com-Youtube.html) | YT | 32 | 1,116,052 | 2,894,183 | Constructed user partitions | Star |
| [AMiner](https://ericdongyx.github.io/metapath2vec/m2v.html) | AM | 3 | 4,891,819 | 12,518,010 | Author, Paper, Venue | Star |

WordNet uses version 3.0, LastFM uses HetRec 2011 tagging events, and Hetionet uses version 1.0. DRKG is a biomedical knowledge graph that includes Hetionet data. BibSonomy connects users, tags, and bookmarked items; its processed input comes from [KONECT](http://konect.cc/). Freebase uses six music entity types from Freebase Easy. YouTube comes from SNAP and uses constructed partitions, not native semantic types. AMiner comes from the metapath2vec dataset release. The other six processed graphs are not bundled.

## Included inputs

WordNet is provided as the full processed 25-layer graph and the three-layer workload used in the main experiments. LastFM is provided as the complete processed three-layer graph.

| File | Contents |
| --- | --- |
| `wordnet_k25.txt` | Relations across 25 selected WordNet lexname classes. |
| `wordnet-p3.txt` | Complete projection onto original layers 12, 13, and 14, renumbered 0, 1, and 2: 2,696 vertices and 1,685 edges, including 21 non-pattern edges. |
| `wordnet-p3.vertices.tsv` | Explicit WordNet-P3 vertex domain; required to retain vertices omitted by the edge list. |
| `lastfm_hetrec_k3.txt` | Pairwise projection of `user_taggedartists.dat` onto User--Artist, User--Tag, and Artist--Tag edges. |

LastFM layers 0, 1, and 2 contain 1,892 users, 12,523 artists, and 9,749 tags. User--Artist edges describe tagging events, not listening counts. WordNet-P3 is a complete layer projection, not a sampled witness; its edge list alone contains 1,483 vertices.

Graph rows contain four whitespace-separated fields:

```text
layer_u vertex_id_u layer_v vertex_id_v
```

Vertex-table rows contain `layer vertex_id`. Layers are 0-based and vertex IDs are local to a layer. Graphs are unweighted, undirected, and simple. Blank lines and lines beginning with `#` are ignored.

## Fixed paper workloads

The experiments keep a fixed layer pattern and retain all input edges. Only pattern edges impose completeness constraints and contribute to the edge objective. Construction uses per-layer thresholds of two; the main Global and Local searches use four, while threshold experiments vary the query thresholds.

WordNet-P3 uses path `0--1--2`, which is also a star. Always supply `wordnet-p3.vertices.tsv` with this input. LastFM and BibSonomy use User-centered stars. Hetionet, Freebase, and YouTube use stars centered on Gene, Genre, and layer 0, respectively.

DRKG uses edges `{(0,2), (0,3), (0,4), (1,4)}` in the layer order listed above. AMiner-AP retains the Author--Paper pair, with 4,887,936 vertices and 9,323,605 edges. Figure abbreviations WN, DK, and AM refer to WordNet-P3, the DRKG tree, and AMiner-AP. The BibSonomy counts above follow the paper's 3,754,828-edge input, not the separate 7,665,240-edge conversion.

See the [project Readme](../Readme.md) for construction and query commands.
