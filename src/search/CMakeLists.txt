add_library(mcviz_search STATIC
    pattern_bounds.cpp pattern_core.cpp pattern_core_queue.cpp
    pattern_solver.cpp pattern_cover.cpp pattern_rules.cpp pattern_pivot.cpp
    star_solver.cpp star_build.cpp star_support_build.cpp
    core_capacity.cpp pattern_filter.cpp
    pattern_dynamic_build.cpp pattern_adjacency.cpp pattern_edge_domains.cpp pattern_partner.cpp)
target_link_libraries(mcviz_search PUBLIC mcviz_common)
target_link_libraries(mcviz_search PRIVATE Threads::Threads)
target_compile_features(mcviz_search PUBLIC cxx_std_17)
