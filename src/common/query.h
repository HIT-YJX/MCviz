#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace mcviz {
using Layer = unsigned;
using VertexId = std::uint64_t;
using VertexPos = unsigned;
using Score = std::uint64_t;
using Thresholds = std::vector<Score>;
using Sets = std::vector<std::vector<VertexPos>>;

struct VertexRef { Layer layer = 0; VertexId id = 0; };
struct LocalVertex { Layer layer = 0; VertexPos pos = 0; };
using Required = std::vector<LocalVertex>;

struct QueryRequest {
    Thresholds gamma;
    std::optional<VertexRef> required;
};

enum class Status { optimal, infeasible, incomplete, error, unsupported };
struct SearchResult {
    bool found = false;
    Score score = 0;
    Sets sets;
    Status status = Status::incomplete;
    std::string message;
};

const char* status_name(Status status);
}
