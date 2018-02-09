#ifndef __KMER_CONTEXT_HPP__
#define __KMER_CONTEXT_HPP__

#include "SalmonUtils.hpp"
#include "UtilityFunctions.hpp"
#include <iostream>
#include <limits>

class KmerContext {
private:
  uint32_t _invalidIdx{std::numeric_limits<uint32_t>::max()};
public:
  KmerContext(uint32_t K, salmon::utils::Direction dir)
      : _valid(false), _length(K), _dir(dir), _repr(_invalidIdx) {
    if (K == 0) {
      std::cerr << "Cannot create a k-mer context of size 0!\n";
      std::exit(1);
    }
  }
  bool valid() const { return _valid; }
  uint32_t index() const { return _repr; }
  std::string str() const {
    return (valid() ? kmerForIndex(_repr, _length) : std::string(_length, 'X'));
  }

  void operator()(const char* s) {
    if (valid()) {
      _repr = nextKmerIndex(_repr, s[_length - 1], _length, _dir);
    } else {
      _repr = 0;
      _repr = indexForKmer(s, _length, _dir);
      _valid = true;
    }
  }

private:
  bool _valid;
  uint32_t _length;
  salmon::utils::Direction _dir;
  uint32_t _repr;
};

#endif //__KMER_CONTEXT_HPP__
