/*  This file is part of Jellyfish.

    Jellyfish is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Jellyfish is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Jellyfish.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __JELLYFISH_BLOOM_COMMON_HPP__
#define __JELLYFISH_BLOOM_COMMON_HPP__

#include <math.h>
#include <jellyfish/divisor.hpp>
#include <jellyfish/atomic_gcc.hpp>

namespace jellyfish {
template<typename Key>
struct hash_pair { };

template<typename Key, typename Derived, typename HashPair = hash_pair<Key> >
class bloom_base {
protected:
  struct prefetch_info {
    size_t         boff;
    unsigned char* pos;
  };

  // The number of bits in the structure, previously known as m_, is
  // know stored as d_.d()
  const jflib::divisor64 d_;
  const unsigned long    k_;
  unsigned char * const  data_;
  HashPair               hash_fns_;

public:
  typedef Key key_type;

  bloom_base(size_t m, unsigned long k, unsigned char* ptr, const HashPair& fns = HashPair()) :
    d_(m), k_(k), data_(ptr), hash_fns_(fns)
  { }

  bloom_base(const bloom_base& rhs) = delete;
  bloom_base(bloom_base&& rhs) :
    d_(rhs.d_), k_(rhs.k_), data_(rhs.data_), hash_fns_(std::move(rhs.hash_fns_))
  { }


  void write_bits(std::ostream& out) {
    out.write((char*)data_, static_cast<Derived*>(this)->nb_bytes());
  }

  // Number of hash functions
  unsigned long k() const { return k_; }
  // Size of bit vector
  size_t m() const { return d_.d(); }
  const HashPair& hash_functions() const { return hash_fns_; }

  static const double LOG2;
  static const double LOG2_SQ;

  static size_t opt_m(const double fp, const size_t n) {
    return n * (size_t)lrint(-log(fp) / LOG2_SQ);
  }
  static unsigned long opt_k(const double fp) {
    return lrint(-log(fp) / LOG2);
  }

  // Insert key k. Returns previous value of k
  unsigned int insert(const Key &k) {
    uint64_t hashes[2];
    hash_fns_(k, hashes);
    return static_cast<Derived*>(this)->insert__(hashes);
  }

  unsigned int check(const Key &k) const {
    uint64_t hashes[2];
    hash_fns_(k, hashes);
    return static_cast<const Derived*>(this)->check__(hashes);
  }



  // Limited std::map interface compatibility
  class element_proxy {
    Derived&   bc_;
    const Key& k_;

  public:
    element_proxy(Derived& bc, const Key& k) : bc_(bc), k_(k) { }

    unsigned int operator++() {
      unsigned int res = bc_.insert(k_);
      return res == 0 ? 1 : 2;
    }

    unsigned int operator++(int) { return bc_.insert(k_); }
    unsigned int operator*() const { return bc_.check(k_); }
    operator unsigned int() const { return bc_.check(k_); }
  };

  class const_element_proxy {
    const Derived& bc_;
    const Key&     k_;

  public:
    const_element_proxy(const Derived& bc, const Key& k) : bc_(bc), k_(k) { }

    unsigned int operator*() const { return bc_.check(k_); }
    operator unsigned int() const { return bc_.check(k_); }
  };
  element_proxy operator[](const Key& k) { return element_proxy(*static_cast<Derived*>(this), k); }
  const_element_proxy operator[](const Key& k) const { return const_element_proxy(*static_cast<const Derived*>(this), k); }
};
template<typename Key, typename Derived, typename HashPair>
const double bloom_base<Key, Derived, HashPair>::LOG2 = 0.6931471805599453;
template<typename Key, typename Derived, typename HashPair>
const double bloom_base<Key, Derived, HashPair>::LOG2_SQ = 0.4804530139182014;

}

#endif /* __JELLYFISH_BLOOM_COMMON_HPP__ */
