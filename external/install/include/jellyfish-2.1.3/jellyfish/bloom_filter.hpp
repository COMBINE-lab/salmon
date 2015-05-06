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

#ifndef __JELLYFISH_BLOOM_FILTER_HPP__
#define __JELLYFISH_BLOOM_FILTER_HPP__

#include <jellyfish/allocators_mmap.hpp>
#include <jellyfish/atomic_gcc.hpp>
#include <jellyfish/bloom_common.hpp>

namespace jellyfish {
template<typename Key, typename HashPair = hash_pair<Key>, typename atomic_t = ::atomic::gcc>
class bloom_filter_base :
    public bloom_base<Key, bloom_filter_base<Key, HashPair>, HashPair>
{
  typedef bloom_base<Key, bloom_filter_base<Key, HashPair>, HashPair> super;

protected:
  static size_t nb_bytes__(size_t l) {
    return l / 8 + (l % 8 != 0);
  }

public:
  bloom_filter_base(size_t m, unsigned long k, unsigned char* ptr, const HashPair& fns = HashPair()) :
    super(m, k, ptr, fns)
  { }

  bloom_filter_base(bloom_filter_base&& rhs) :
  super(std::move(rhs))
  { }

  size_t nb_bytes() const {
    return nb_bytes__(super::d_.d());
  }

  // Insert key with given hashes
  unsigned int insert__(const uint64_t *hashes) {
    // Prefetch memory locations
    // This static_assert make clang++ happy...
    static_assert(std::is_pod<typename super::prefetch_info>::value, "prefetch_info must be a POD");

    typename super::prefetch_info pinfo[super::k_];
    const size_t base    = super::d_.remainder(hashes[0]);
    const size_t inc     = super::d_.remainder(hashes[1]);
    for(unsigned long i = 0; i < super::k_; ++i) {
      const size_t pos   = super::d_.remainder(base + i * inc);
      const size_t elt_i = pos / 8;
      pinfo[i].boff      = pos % 8;
      pinfo[i].pos       = super::data_ + elt_i;
      __builtin_prefetch(pinfo[i].pos, 1, 0);
    }

    // Check if element present
    bool present = true;
    for(unsigned long i = 0; i < super::k_; ++i) {
      const char mask = (char)1 << pinfo[i].boff;
      const char prev = __sync_fetch_and_or(pinfo[i].pos, mask);
      present         = present && (prev & mask);
    }

    return present;
  }

  // Compute hashes of key k
  void hash(const Key &k, uint64_t *hashes) const { hash_fns_(k, hashes); }

  unsigned int check__(const uint64_t *hashes) const {
    // Prefetch memory locations
    typename super::prefetch_info pinfo[super::k_];
    const size_t base    = super::d_.remainder(hashes[0]);
    const size_t inc     = super::d_.remainder(hashes[1]);
    for(unsigned long i = 0; i < super::k_; ++i) {
      const size_t pos   = super::d_.remainder(base + i * inc);
      const size_t elt_i = pos / 8;
      pinfo[i].boff      = pos % 8;
      pinfo[i].pos       = super::data_ + elt_i;
      __builtin_prefetch(pinfo[i].pos, 0, 0);
    }

    for(unsigned long i = 0; i < super::k_; ++i)
      if(!(jflib::a_load(pinfo[i].pos) & ((char)1 << pinfo[i].boff)))
        return 0;
    return 1;
  }
};

template<typename Key, typename HashPair = hash_pair<Key>, typename atomic_t = ::atomic::gcc,
         typename mem_block_t = allocators::mmap>
class bloom_filter :
  protected mem_block_t,
  public bloom_filter_base<Key, HashPair, atomic_t>
{
  typedef bloom_filter_base<Key, HashPair, atomic_t> super;
public:
  bloom_filter(double fp, size_t n, const HashPair& fns = HashPair()) :
    mem_block_t(super::nb_bytes__(super::opt_m(fp, n))),
    super(super::opt_m(fp, n), super::opt_k(fp), (unsigned char*)mem_block_t::get_ptr(), fns)
  {
    if(!mem_block_t::get_ptr())
      eraise(std::runtime_error) << "Failed to allocate " << super::nb_bytes__(super::opt_m(fp, n))
                                 << " bytes of memory for bloom_filter";
  }

  bloom_filter(size_t m, unsigned long k, const HashPair& fns = HashPair()) :
    mem_block_t(super::nb_bytes__(m)),
    super(m, k, (unsigned char*)mem_block_t::get_ptr(), fns)
  {
    if(!mem_block_t::get_ptr())
      eraise(std::runtime_error) << "Failed to allocate " << super::nb_bytes__(m) << " bytes of memory for bloom_filter";
  }

  bloom_filter(size_t m, unsigned long k, std::istream& is, const HashPair& fns = HashPair()) :
    mem_block_t(super::nb_bytes__(m)),
    super(m, k, (unsigned char*)mem_block_t::get_ptr(), fns)
  {
    if(!mem_block_t::get_ptr())
      eraise(std::runtime_error) << "Failed to allocate " << super::nb_bytes__(m) << " bytes of memory for bloom_filter";

    is.read((char*)mem_block_t::get_ptr(), mem_block_t::get_size());
  }

  bloom_filter(bloom_filter&& rhs) :
    mem_block_t(std::move(rhs)),
    super(std::move(rhs))
  { }
};

template<typename Key, typename HashPair = hash_pair<Key>, typename atomic_t = ::atomic::gcc>
class bloom_filter_file :
    protected mapped_file,
    public bloom_filter_base<Key, HashPair, atomic_t>
{
  typedef bloom_filter_base<Key, HashPair, atomic_t> super;
public:
  typedef typename super::key_type key_type;

  bloom_filter_file(size_t m, unsigned long k, const char* path, const HashPair& fns = HashPair(), off_t offset = 0) :
    mapped_file(path),
    super(m, k, (unsigned char*)mapped_file::base() + offset, fns)
  { }

  bloom_filter_file(const bloom_filter_file& rhs) = delete;
  bloom_filter_file(bloom_filter_file&& rhs) :
    mapped_file(std::move(rhs)),
    super(std::move(rhs))
  { }
};

}

#endif /* __JELLYFISH_BLOOM_FILTER_HPP__ */
