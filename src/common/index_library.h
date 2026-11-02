#pragma once
#include "cli.h"

namespace mcviz {
struct LibraryPaths {
    std::string index, stem, prefix;
};
LibraryPaths library_paths(const CliOptions& options, const Graph& graph,
                           const BuildSpec& spec, const std::string& algorithm);


BuildIdentity library_index_identity(const std::string& path);
std::string find_library_index(const LibraryPaths& paths, const Graph& graph,
                               const BuildSpec& spec, const std::string& engine,
                               const std::string& value_kind = {});
}
