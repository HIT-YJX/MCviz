#include "text.h"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace mcviz {

void require(bool ok, const std::string& message) {
    if (!ok) throw std::runtime_error(message);
}


std::uint64_t number(const std::string& s, const std::string&, std::uint64_t) {
    std::uint64_t n = 0;
    for (char c : s) n = n * 10 + static_cast<unsigned>(c - '0');
    return n;
}


std::string trim(const std::string& s) {
    auto a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return {};
    return s.substr(a, s.find_last_not_of(" \t\r\n") - a + 1);
}


std::string read_bytes(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    std::ostringstream out;
    out << in.rdbuf();
    return out.str();
}


std::vector<std::string> data_lines(const std::string& bytes) {
    std::istringstream in(bytes);
    std::vector<std::string> rows;
    std::string s;
    while (std::getline(in, s)) {
        s = trim(s.substr(0, s.find('#')));
        if (!s.empty()) rows.push_back(s);
    }
    return rows;
}


std::vector<std::string> fields(const std::string& line) {
    std::istringstream in(line);
    std::vector<std::string> out;
    std::string s;
    while (in >> s) out.push_back(s);
    return out;
}


std::vector<std::string> csv_fields(const std::string& line) {
    std::vector<std::string> out;
    std::size_t start = 0;
    while (true) {
        auto end = line.find(',', start);
        out.push_back(trim(line.substr(start, end == std::string::npos ? end : end - start)));
        if (end == std::string::npos) return out;
        start = end + 1;
    }
}


const char* status_name(Status s) {
    switch (s) {
    case Status::optimal: return "optimal";
    case Status::infeasible: return "infeasible";
    case Status::incomplete: return "incomplete";
    case Status::error: return "error";
    case Status::unsupported: return "unsupported";
    }
    throw std::runtime_error("Invalid result status");
}
}
