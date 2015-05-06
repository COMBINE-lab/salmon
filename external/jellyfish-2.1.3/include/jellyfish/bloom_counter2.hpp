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


#ifndef __BLOOM_COUNTER2_HPP__
#define __BLOOM_COUNTER2_HPP__

#include <assert.h>
#include <math.h>
#include <limits.h>
#include <jellyfish/bloom_common.hpp>
#include <jellyfish/mapped_file.hpp>
#include <jellyfish/allocators_mmap.hpp>
#include <jellyfish/divisor.hpp>
#include <jellyfish/atomic_gcc.hpp>
#include <jellyfish/atomic_field.hpp>

#include <jellyfish/err.hpp>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

namespace jellyfish {
/* Bloom counter with 3 values: 0, 1 or 2. It is thread safe and lock free.
 */
template<typename Key, typename HashPair = hash_pair<Key>, typename atomic_t = ::atomic::gcc>
class bloom_counter2_base : public bloom_base<Key, bloom_counter2_base<Key, HashPair, atomic_t>, HashPair> {
  typedef bloom_base<Key, bloom_counter2_base<Key, HashPair, atomic_t>, HashPair> super;

  atomic_t atomic_;


protected:
  static size_t nb_bytes__(size_t l) {
    return l / 5 + (l % 5 != 0);
  }

public:
  bloom_counter2_base(size_t m, unsigned long k, unsigned char* ptr, const HashPair& fns = HashPair()) :
    super(m, k, ptr, fns)
  { }
  bloom_counter2_base(bloom_counter2_base&& rhs) :
    super(std::move(rhs))
  { }
  size_t nb_bytes() const {
    return nb_bytes__(super::d_.d());
  }

  // Insert key with given hashes
  unsigned int insert__(const uint64_t* hashes) {
    // Prefetch memory locations
    static_assert(std::is_pod<typename super::prefetch_info>::value, "prefetch_info must be a POD");
    typename super::prefetch_info pinfo[super::k_];
    const size_t base = super::d_.remainder(hashes[0]);
    const size_t inc  = super::d_.remainder(hashes[1]);
    for(unsigned long i = 0; i < super::k_; ++i) {
      const size_t p   = super::d_.remainder(base + i * inc);
      const size_t off = p / 5;
      pinfo[i].boff    = p % 5;
      pinfo[i].pos     = super::data_ + off;
      //      prefetch_write_no(pinfo[i].pos);
      __builtin_prefetch(pinfo[i].pos, 1, 0);
    }

    // Insert element
    unsigned char res = 2;
    for(unsigned long i = 0; i < super::k_; ++i) {
      size_t        boff = pinfo[i].boff;
      unsigned char v    = jflib::a_load(pinfo[i].pos);

      while(true) {
        unsigned char w = v;
        switch(boff) {
        case 0:          break;
        case 1: w /= 3;  break;
        case 2: w /= 9;  break;
        case 3: w /= 27; break;
        case 4: w /= 81; break;
        }
        w = w % 3;
        if(w == 2) break;
        unsigned char nv = v;

        switch(boff) {
        case 0: nv += 1;  break;
        case 1: nv += 3;  break;
        case 2: nv += 9;  break;
        case 3: nv += 27; break;
        case 4: nv += 81; break;
        }
        unsigned char cv = atomic_.cas(pinfo[i].pos, v, nv);
        if(cv == v) {
          if(w < res)
            res = w;
          break;
        }
        v = cv;
      }
    }
    return res;
  }

  unsigned int check__(uint64_t *hashes) const {
    // Prefetch memory locations
    static_assert(std::is_pod<typename super::prefetch_info>::value, "prefetch_info must be a POD");
    typename super::prefetch_info pinfo[super::k_];
    const size_t base = super::d_.remainder(hashes[0]);
    const size_t inc  = super::d_.remainder(hashes[1]);
    for(unsigned long i = 0; i < super::k_; ++i) {
      const size_t p   = super::d_.remainder(base + i * inc);
      const size_t off = p / 5;
      pinfo[i].boff    = p % 5;
      pinfo[i].pos     = super::data_ + off;
      //      prefetch_read_no(pinfo[i].pos);
      __builtin_prefetch(pinfo[i].pos, 0, 0);
    }

    // Check element
    unsigned char res = 2;
    for(unsigned long i = 0; i < super::k_; ++i) {
      size_t        boff = pinfo[i].boff;
      unsigned char w    = jflib::a_load(pinfo[i].pos);

      switch(boff) {
      case 0:          break;
      case 1: w /= 3;  break;
      case 2: w /= 9;  break;
      case 3: w /= 27; break;
      case 4: w /= 81; break;
      }
      w = w % 3;
      if(w < res)
        res = w;
    }
    return res;
  }
};

template<typename Key, typename HashPair = hash_pair<Key>, typename atomic_t = ::atomic::gcc,
         typename mem_block_t = allocators::mmap>
class bloom_counter2:
    protected mem_block_t,
    public bloom_counter2_base<Key, HashPair, atomic_t>
{
  typedef bloom_counter2_base<Key, HashPair, atomic_t> super;

public:
  typedef typename super::key_type key_type;

  bloom_counter2(const double fp, const size_t n, const HashPair& fns = HashPair()) :
    mem_block_t(super::nb_bytes__(super::opt_m(fp, n))),
    super(super::opt_m(fp, n), super::opt_k(fp), (unsigned char*)mem_block_t::get_ptr(), fns)
  {
    if(!mem_block_t::get_ptr())
      eraise(std::runtime_error) << "Failed to allocate " << super::nb_bytes__(super::opt_m(fp, n))
                                 << " bytes of memory for bloom_counter";
  }

  bloom_counter2(size_t m, unsigned long k, const HashPair& fns = HashPair()) :
    mem_block_t(super::nb_bytes__(m)),
    super(m, k, (unsigned char*)mem_block_t::get_ptr(), fns)
  {
    if(!mem_block_t::get_ptr())
      eraise(std::runtime_error) << "Failed to allocate " << super::nb_bytes__(m) << " bytes of memory for bloom_counter";
  }

  bloom_counter2(size_t m, unsigned long k, std::istream& is, const HashPair& fns = HashPair()) :
    mem_block_t(super::nb_bytes__(m)),
    super(m, k, (unsigned char*)mem_block_t::get_ptr(), fns)
  {
    if(!mem_block_t::get_ptr())
      eraise(std::runtime_error) << "Failed to allocate " << super::nb_bytes__(m) << " bytes of memory for bloom_counter";

    is.read((char*)mem_block_t::get_ptr(), mem_block_t::get_size());
  }

  bloom_counter2(const bloom_counter2& rhs) = delete;
  bloom_counter2(bloom_counter2&& rhs) :
    mem_block_t(std::move(rhs)),
    super(std::move(rhs))
  { }
};

template<typename Key, typename HashPair = hash_pair<Key>, typename atomic_t = ::atomic::gcc>
class bloom_counter2_file :
    protected mapped_file,
    public bloom_counter2_base<Key, HashPair, atomic_t>
{
  typedef bloom_counter2_base<Key, HashPair, atomic_t> super;
public:
  typedef typename super::key_type key_type;

  bloom_counter2_file(size_t m, unsigned long k, const char* path, const HashPair& fns = HashPair(), off_t offset = 0) :
    mapped_file(path),
    super(m, k, (unsigned char*)mapped_file::base() + offset, fns)
  { }

  bloom_counter2_file(const bloom_counter2_file& rhs) = delete;
  bloom_counter2_file(bloom_counter2_file&& rhs) :
    mapped_file(std::move(rhs)),
    super(std::move(rhs))
  { }
};

} // namespace jellyfish {

#endif // __BLOOM_COUNTER2_HPP__
