#pragma once
#include "cli.h"
#include <cstdint>

namespace mcviz {
void log_build(const CliOptions& opt, const BuildSpec& spec, const std::string& method,
               double seconds, std::uintmax_t bytes, bool reused = false);
void log_search(const CliOptions& opt, const BuildSpec& spec, const QueryRequest& query,
                const SearchResult& result, const std::string& method, double seconds);
void log_conversion(const CliOptions& opt, const BuildSpec& spec, const std::string& method,
                    double forward, double reverse, std::uint64_t tree_bytes,
                    std::uint64_t linear_bytes);
}
