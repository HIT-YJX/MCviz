#ifndef BASELINE_COMMON_H
#define BASELINE_COMMON_H

#include "../common/graph.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <chrono>
#include <cmath>

using namespace std;
using namespace std::chrono;

typedef unsigned int uint;
typedef unsigned long long ull;

enum ProblemType {
    MESC = 1,
    MVSC = 2,
    MBSC = 3
};

struct StarClique {
    vector<vector<uint>> sets;
    ull size;
    bool found;

    StarClique() : size(0), found(false) {}
    StarClique(uint K) : sets(K), size(0), found(false) {}
};

#endif
