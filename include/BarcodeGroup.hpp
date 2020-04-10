#ifndef BARCODE_GROUP_HPP
#define BARCODE_GROUP_HPP

#include <cuckoohash_map.hh>
#include "xxhash.h"
#include "metro/metrohash.h"
#include "tsl/array_map.h"

struct BarcodeGroupStringHasher {
  std::size_t operator()(const std::string& k) const {
    char* pt = const_cast<char*>(k.data());
    return XXH64(static_cast<void*>(pt), k.size() * sizeof(char), 0);
  }
};

struct BarcodeGroupStringHasherMetro {
  std::size_t operator()(const char* key, std::size_t key_size) const {
    std::size_t r;
    MetroHash64::Hash(reinterpret_cast<const uint8_t*>(key), key_size, reinterpret_cast<uint8_t*>(&r), 0);
    return r;
  }
};

//using CFreqMapT = libcuckoo::cuckoohash_map<std::string, uint32_t, BarcodeGroupStringHasher>;
using CFreqMapT = tsl::array_map<char, uint32_t, BarcodeGroupStringHasherMetro>;

using MapT = std::vector<std::pair<std::string, double>>;
using TrueBcsT = std::unordered_set<std::string>;
using TrueMapT = std::unordered_map<std::string, uint32_t>;
using SoftMapT = std::unordered_map<std::string, MapT>;

#endif // BARCODE_GROUP_HPP
