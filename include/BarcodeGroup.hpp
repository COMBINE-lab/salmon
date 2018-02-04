#ifndef BARCODE_GROUP_HPP
#define BARCODE_GROUP_HPP

#include <cuckoohash_map.hh>
#include "xxhash.h"

struct BarcodeGroupStringHasher {
  std::size_t operator()(const std::string& k) const {
    char* pt = const_cast<char*>(k.data());
    return XXH64(static_cast<void*>(pt), k.size() * sizeof(char), 0);
  }
};

using CFreqMapT = cuckoohash_map<std::string, uint32_t, BarcodeGroupStringHasher>;

using MapT = std::vector<std::pair<std::string, double>>;
using TrueBcsT = std::unordered_set<std::string>;
using TrueMapT = std::unordered_map<std::string, uint32_t>;
using SoftMapT = std::unordered_map<std::string, MapT>;

#endif // BARCODE_GROUP_HPP
