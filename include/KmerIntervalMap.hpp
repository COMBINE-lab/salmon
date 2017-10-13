#ifndef __KMER_INTERVAL_MAP_HPP__
#define __KMER_INTERVAL_MAP_HPP__

extern "C" {
#include "bwt.h"
}

#include <fstream>
#include <unordered_map>

#include <jellyfish/mer_dna.hpp>

#include "cereal/archives/binary.hpp"
#include "cereal/types/unordered_map.hpp"

#include "xxhash.h"

using JFMer = jellyfish::mer_dna_ns::mer_base_static<uint64_t, 1>;

// What will be the keys in our k-mer has map
struct KmerKey {
  KmerKey() { mer_.polyT(); }

  KmerKey(uint8_t* seq, uint32_t len) : mer_(len) {
    mer_.polyT();
    for (size_t i = 0; i < len; ++i) {
      mer_.shift_left(seq[i]);
    }
  }

  bool operator==(const KmerKey& ok) const { return mer_ == ok.mer_; }

  // Is there a smarter way to do save / load here?
  template <typename Archive> void save(Archive& archive) const {
    auto key = mer_.get_bits(0, 2 * mer_.k());
    archive(key);
  }

  template <typename Archive> void load(Archive& archive) {
    mer_.polyT();
    uint64_t bits;
    archive(bits);
    mer_.set_bits(0, 2 * mer_.k(), bits);
  }

  JFMer mer_;
};

template <typename Archive> void load(Archive& archive, bwtintv_t& interval) {
  archive(interval.x[0], interval.x[1], interval.x[2], interval.info);
}

template <typename Archive>
void save(Archive& archive, const bwtintv_t& interval) {
  archive(interval.x[0], interval.x[1], interval.x[2], interval.info);
}

/**
 *  This class provides an efficent hash-map from
 *  k-mers to BWT intervals.
 */
class KmerIntervalMap {
public:
  // How we hash the keys
  struct KmerHasher {
    std::size_t operator()(const KmerKey& k) const {
      void* data = static_cast<void*>(const_cast<KmerKey&>(k).mer_.data__());
      return XXH64(data, sizeof(uint64_t), 0);
    }
  };

private:
  std::unordered_map<KmerKey, bwtintv_t, KmerHasher> map_;

public:
  void setK(unsigned int k) { JFMer::k(k); }
  uint32_t k() { return JFMer::k(); }

  bool hasKmer(KmerKey& k) { return map_.find(k) != map_.end(); }

  decltype(map_)::iterator find(const KmerKey& k) { return map_.find(k); }
  decltype(map_)::iterator find(KmerKey&& k) { return map_.find(k); }

  decltype(map_)::iterator end() { return map_.end(); }

  bwtintv_t& operator[](const KmerKey& k) { return map_[k]; }
  bwtintv_t& operator[](KmerKey&& k) { return map_[k]; }

  decltype(map_)::size_type size() { return map_.size(); }

  void save(boost::filesystem::path indexPath) {
    std::ofstream ofs(indexPath.string(), std::ios::binary);
    {
      cereal::BinaryOutputArchive oa(ofs);
      oa(map_);
    }
    ofs.close();
  }

  void load(boost::filesystem::path indexPath) {
    std::ifstream ifs(indexPath.string(), std::ios::binary);
    {
      cereal::BinaryInputArchive ia(ifs);
      ia(map_);
    }
    ifs.close();
  }
};

#endif // __KMER_INTERVAL_MAP_HPP__
