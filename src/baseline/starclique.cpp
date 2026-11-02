#include "starclique.h"

bool StarCliqueFinder::contains_required(const vector<vector<uint>>& sets) const {
    if (!required) return true;
    const auto& layer = sets[required->layer];
    return find(layer.begin(), layer.end(), required->pos) != layer.end();
}
