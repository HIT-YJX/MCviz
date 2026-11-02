add_library(mcviz_common STATIC
    text.cpp hash.cpp graph.cpp pattern.cpp objective.cpp contract.cpp cli.cpp log.cpp index_library.cpp)
target_include_directories(mcviz_common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_compile_features(mcviz_common PUBLIC cxx_std_17)
target_compile_definitions(mcviz_common PUBLIC MCVIZ_RELEASE_ROOT="${CMAKE_SOURCE_DIR}")
