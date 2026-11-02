add_library(mcviz_pattern_index STATIC
    forest.cpp pattern_index.cpp pattern_query.cpp query_seed_core.cpp query_seed_heuristic.cpp
    state.cpp engine_cli.cpp build_filter.cpp)
target_link_libraries(mcviz_pattern_index PUBLIC mcviz_search)
target_compile_features(mcviz_pattern_index PUBLIC cxx_std_17)
add_executable(mcviz_representation_conversion
    conversion_main.cpp representation.cpp)
target_link_libraries(mcviz_representation_conversion PRIVATE mcviz_pattern_index)
