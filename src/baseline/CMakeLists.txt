add_executable(baseline_layer
    main.cpp api.cpp bcvizm_lp.cpp bcvizm_order.cpp bcvizm_index.cpp kgraph_adapter.cpp
    starclique.cpp starclique_motif.cpp
    starclique_kpc.cpp starclique_kspc.cpp starclique_kspc_star.cpp)
target_link_libraries(baseline_layer PRIVATE mcviz_common)
