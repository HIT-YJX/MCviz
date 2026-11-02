#pragma once

#include "query.h"
#include <limits>
#include <string>
#include <vector>

namespace mcviz {
void require(bool ok, const std::string& message);
std::uint64_t number(const std::string& s, const std::string& field,
                     std::uint64_t limit = std::numeric_limits<std::uint64_t>::max());
std::string trim(const std::string& s);
std::string read_bytes(const std::string& path);
std::vector<std::string> data_lines(const std::string& bytes);
std::vector<std::string> fields(const std::string& line);
std::vector<std::string> csv_fields(const std::string& line);
std::string sha256(const std::string& bytes);
}
