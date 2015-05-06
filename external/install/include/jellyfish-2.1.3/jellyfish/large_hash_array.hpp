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

#ifndef __JELLYFISH_LARGE_HASH_ARRAY_HPP__
#define __JELLYFISH_LARGE_HASH_ARRAY_HPP__

#include <jellyfish/storage.hpp>
#include <jellyfish/atomic_gcc.hpp>
#include <jellyfish/allocators_mmap.hpp>
#include <jellyfish/offsets_key_value.hpp>
#include <jellyfish/misc.hpp>
#include <jellyfish/err.hpp>
#include <jellyfish/mer_dna.hpp>
#include <jellyfish/rectangular_binary_matrix.hpp>
#include <jellyfish/simple_circular_buffer.hpp>
#include <jellyfish/large_hash_iterator.hpp>

namespace jellyfish { namespace large_hash {
/* Contains an integer, the reprobe limit. It is capped based on the
 * reprobe strategy to not be bigger than the size of the hash
 * array. Also, the length to encode reprobe limit must not be larger
 * than the length to encode _size.
 */
class reprobe_limit_t {
  uint_t limit;
public:
  reprobe_limit_t(uint_t _limit, const size_t *_reprobes, size_t _size) :
    limit(_limit)
  {
    while(_reprobes[limit] >= _size && limit >= 1)
      limit--;
  }
  inline uint_t val() const { return limit; }
};

// Key is any type with the following two methods: get_bits(unsigned
// int start, unsigned int len); and set_bits(unsigned int start,
// unsigned int len, uint64_t bits). These methods get and set the
// bits [start, start + len). Start and len may not be aligned to word
// boundaries. On the other hand, len is guaranteed to be <
// sizeof(uint64_t). I.e. never more than 1 word is fetched or set.
template<typename Key, typename word, typename atomic_t, typename Derived>
class array_base {
  static const int  wsize = std::numeric_limits<word>::digits; // Word size in bits
  // Can't be done. Resort to an evil macro!
  //  static const word fmask = std::numeric_limits<word>::max(); // Mask full of ones
#define fmask (std::numeric_limits<word>::max())

public:
  define_error_class(ErrorAllocation);

  typedef word                             data_word;
  typedef typename Offsets<word>::offset_t offset_t;
  typedef struct offset_t::key             key_offsets;
  typedef struct offset_t::val             val_offsets;

  typedef Key key_type;
  typedef uint64_t                      mapped_type;
  typedef std::pair<Key&, mapped_type>  value_type;
  typedef stl_iterator_base<array_base> iterator;
  typedef stl_iterator_base<array_base> const_iterator;
  typedef value_type&                   reference;
  typedef const value_type&             const_reference;
  typedef value_type*                   pointer;
  typedef const value_type*             const_pointer;

  typedef eager_iterator_base<array_base>  eager_iterator;
  typedef lazy_iterator_base<array_base>   lazy_iterator;
  typedef region_iterator_base<array_base> region_iterator;

  /// Status of a (key,value) pair. LBSET means that the large bit is
  /// set. Hence, it contains a pointer back to the original key and a
  /// large value.
  enum key_status { FILLED, EMPTY, LBSET};

protected:
  uint16_t                 lsize_; // log of size
  size_t                   size_, size_mask_;
  reprobe_limit_t          reprobe_limit_;
  uint16_t                 key_len_; // Length of key in bits
  uint16_t                 raw_key_len_; // Length of key stored raw (i.e. complement of implied length)
  Offsets<word>            offsets_; // key len reduced by size of hash array
  size_t                   size_bytes_;
  word * const             data_;
  atomic_t                 atomic_;
  const size_t            *reprobes_;
  RectangularBinaryMatrix  hash_matrix_;
  RectangularBinaryMatrix  hash_inverse_matrix_;

public:
  /// Give information about memory usage and array size.
  struct usage_info {
    uint16_t      key_len_, val_len_, reprobe_limit_;
    const size_t* reprobes_;

    usage_info(uint16_t key_len, uint16_t val_len, uint16_t reprobe_limit,
               const size_t* reprobes = jellyfish::quadratic_reprobes) :
      key_len_(key_len), val_len_(val_len), reprobe_limit_(reprobe_limit), reprobes_(reprobes) { }

    /// Memory usage for a given size.
    size_t mem(size_t size) {
      uint16_t lsize(ceilLog2(size));
      size_t asize((size_t)1 << lsize);
      reprobe_limit_t areprobe_limit(reprobe_limit_, reprobes_, asize);
      uint16_t raw_key_len(key_len_ > lsize ? key_len_ - lsize : 0);
      Offsets<word> offsets(raw_key_len + bitsize(areprobe_limit.val() + 1), val_len_,
                            areprobe_limit.val() + 1);
      return div_ceil(asize,
                      (size_t)offsets.block_len()) * offsets.block_word_len() * sizeof(word) + sizeof(array_base) + sizeof(Offsets<word>);
    }

    /// Actual size for a given size.
    size_t asize(size_t size) { return (size_t)1 << ceilLog2(size); }

    struct fit_in {
      usage_info* i_;
      size_t      mem_;
      fit_in(usage_info* i, size_t mem) : i_(i), mem_(mem) { }
      bool operator()(uint16_t size_bits) const { return i_->mem((size_t)1 << size_bits) < mem_; }
    };

    /// Maximum size for a given maximum memory.
    size_t size(size_t mem) { return (size_t)1 << size_bits(mem); }

    /// Log of maximum size for a given maximum memory
    uint16_t size_bits(size_t mem) {
      uint16_t res = *binary_search_first_false(pointer_integer<uint16_t>(0), pointer_integer<uint16_t>(64),
                                              fit_in(this, mem));
      return res > 0 ? res - 1 : 0;
    }

    size_t size_bits_linear(size_t mem) {
      fit_in predicate(this, mem);
      uint16_t i = 0;
      for( ; i < 64; ++i)
        if(!predicate(i))
           break;

      return i > 0 ? i - 1 : 0;
    }

  };


  array_base(size_t size, // Size of hash. To be rounded up to a power of 2
             uint16_t key_len, // Size of key in bits
             uint16_t val_len, // Size of val in bits
             uint16_t reprobe_limit, // Maximum reprobe
             RectangularBinaryMatrix m,
             const size_t* reprobes = quadratic_reprobes) : // Reprobing policy
    lsize_(ceilLog2(size)),
    size_((size_t)1 << lsize_),
    size_mask_(size_ - 1),
    reprobe_limit_(reprobe_limit, reprobes, size_),
    key_len_(key_len),
    raw_key_len_(key_len_ > lsize_ ? key_len_ - lsize_ : 0),
    offsets_(raw_key_len_ + bitsize(reprobe_limit_.val() + 1), val_len, reprobe_limit_.val() + 1),
    size_bytes_(div_ceil(size_, (size_t)offsets_.block_len()) * offsets_.block_word_len() * sizeof(word)),
    data_(static_cast<Derived*>(this)->alloc_data(size_bytes_)),
    reprobes_(reprobes),
    hash_matrix_(m),
    hash_inverse_matrix_(hash_matrix_.pseudo_inverse())
  {
    if(!data_)
      eraise(ErrorAllocation) << "Failed to allocate "
                              << (div_ceil(size, (size_t)offsets_.block_len()) * offsets_.block_word_len() * sizeof(word))
                              << " bytes of memory";
  }

  array_base(array_base&& ary) :
    lsize_(ary.lsize_),
    size_(ary.size_),
    size_mask_(size_ - 1),
    reprobe_limit_(ary.reprobe_limit_),
    key_len_(ary.key_len_),
    raw_key_len_(ary.raw_key_len_),
    offsets_(std::move(ary.offsets_)),
    size_bytes_(ary.size_bytes_),
    data_(ary.data_),
    reprobes_(ary.reprobes_),
    hash_matrix_(std::move(ary.hash_matrix_)),
    hash_inverse_matrix_(std::move(ary.hash_inverse_matrix_))
  { }

  array_base& operator=(const array_base& rhs) = delete;
  array_base& operator=(array_base&& rhs) = delete;

  size_t size() const { return size_; }
  size_t lsize() const { return lsize_; }
  size_t size_mask() const { return size_mask_; }
  uint_t key_len() const { return key_len_; }
  uint_t val_len() const { return offsets_.val_len(); }

  const size_t* reprobes() const { return reprobes_; }
  uint_t max_reprobe() const { return reprobe_limit_.val(); }
  size_t max_reprobe_offset() const { return reprobes_[reprobe_limit_.val()]; }

  const RectangularBinaryMatrix& matrix() const { return hash_matrix_; }
  const RectangularBinaryMatrix& inverse_matrix() const { return hash_inverse_matrix_; }
  void matrix(const RectangularBinaryMatrix& m) {
    hash_inverse_matrix_ = m.pseudo_inverse();
    hash_matrix_         = m;
  }

  /**
   * Clear hash table. Not thread safe.
   */
  void clear() {
    memset(data_, '\0', size_bytes_);
  }

  /**
   * Write the hash table raw to a stream. Not thread safe.
   */
  void write(std::ostream& os) const {
    os.write((const char*)data_, size_bytes_);
  }

  size_t size_bytes() const { return size_bytes_; }

  /* The storage of the hash is organized in "blocks". A (key,value)
   * pair always start at bit 0 of the block. The following methods
   * work with the blocks of the hash.
   */

  /**
   * Number of blocks needed to fit at least a given number of
   * records. Given a number of records, it returns the number of
   * blocks necessary and the actual number of records these blocks
   * contain.
   */
  std::pair<size_t, size_t> blocks_for_records(size_t nb_records) const {
    return offsets_.blocks_for_records(nb_records);
  }


  /**
   * Convert coordinate from (start, blen) given in blocks to
   * coordinate in char* and length in bytes. It also makes sure that
   * the pointer and length returned do not go beyond allocated
   * memory.
   */
  void block_to_ptr(const size_t start, const size_t blen,
                    char **start_ptr, size_t *memlen) const {
    *start_ptr    = (char *)(data_ + start * offsets_.block_word_len());
    char *end_ptr = (char *)data_ + size_bytes_;

    if(*start_ptr >= end_ptr) {
      *memlen = 0;
      return;
    }
    *memlen = blen * offsets_.block_word_len() * sizeof(word);
    if(*start_ptr + *memlen > end_ptr)
      *memlen = end_ptr - *start_ptr;
  }

  /**
   * Zero out blocks in [start, start+length), where start and
   * length are given in number of blocks.
   **/
  void zero_blocks(const size_t start, const size_t length) {
    char   *start_ptr;
    size_t  memlen;
    block_to_ptr(start, length, &start_ptr, &memlen);
    memset(start_ptr, '\0', memlen);
  }


  /**
   * Use hash values as counters.
   *
   * The matrix multiplication gets only a uint64_t. The lsb of the
   * matrix product, the hsb are assume to be equal to the key itself
   * (the matrix has a partial identity on the first rows).
   *
   * In case of failure (false is returned), carry_shift contains the
   * number of bits of the value that were successfully stored in the
   * hash (low significant bits). If carry_shift == 0, then nothing
   * was stored and the key is not in the hash at all. In that case,
   * the value of *is_new and *id are not valid. If carry_shift > 0,
   * then the key is present but the value stored is not correct
   * (missing the high significant bits of value), but *is_new and *id
   * contain the proper information.
   */
  inline bool add(const key_type& key, mapped_type val, unsigned int* carry_shift, bool* is_new, size_t* id) {
    uint64_t hash = hash_matrix_.times(key);
    *carry_shift  = 0;
    return add_rec(hash & size_mask_, key, val, false, is_new, id, carry_shift);
  }

  inline bool add(const key_type& key, mapped_type val, unsigned int* carry_shift) {
    bool   is_new = false;
    size_t id     = 0;
    return add(key, val, carry_shift, &is_new, &id);
  }

  inline bool add(const key_type& key, mapped_type val) {
    unsigned int carry_shift = 0;
    return add(key, val, &carry_shift);
  }

  inline bool set(const key_type& key) {
    bool   is_new;
    size_t id;
    return set(key, &is_new, &id);
  }
  bool set(const key_type& key, bool* is_new, size_t* id) {
    word*           w;
    const offset_t* o;

    *id = hash_matrix_.times(key) & size_mask_;
    return claim_key(key, is_new, id, &o, &w);
  }

  /**
   * Use hash values as counters, if already exists
   *
   * Add val to the value associated with key if key is already in the
   * hash. Returns true if the update was done, false otherwise.
   */
  inline bool update_add(const key_type& key, mapped_type val) {
    key_type     tmp_key;
    unsigned int carry_shift;
    return update_add(key, val, &carry_shift, tmp_key);
  }


  // Optimization. Use tmp_key as buffer. Avoids allocation if update_add is called repeatedly.
  bool update_add(const key_type& key, mapped_type val, unsigned int* carry_shift, key_type& tmp_key) {
    size_t          id;
    word*           w;
    const offset_t* o;
    *carry_shift = 0;

    if(get_key_id(key, &id, tmp_key, (const word**)&w, &o))
      return add_rec_at(id, key, val, o, w, carry_shift);
    return false;
  }

  // Get the value, stored in *val, associated with key. If the key is
  // not found, false is returned, otherwise true is returned and *val
  // is updated. If carry_bit is true, then the first bit of the key
  // field indicates whether we should reprobe to get the complete
  // value.
  inline bool get_val_for_key(const key_type& key, mapped_type* val, bool carry_bit = false) const {
    key_type tmp_key;
    size_t   id;
    return get_val_for_key(key, val, tmp_key, &id, carry_bit);
  }

  // Optimization version. A tmp_key buffer is passed and the id where
  // the key was found is return in *id. If get_val_for_key is called
  // many times consecutively, it may be faster to pass the same
  // tmp_key buffer instead of allocating it every time.
  bool get_val_for_key(const key_type& key, mapped_type* val, key_type& tmp_key,
                       size_t* id, bool carry_bit = false) const {
    const word*     w;
    const offset_t* o;
    if(!get_key_id(key, id, tmp_key, &w, &o))
      return false;
    *val = get_val_at_id(*id, w, o, true, carry_bit);
    return true;
  }

  // Return true if the key is present in the hash
  inline bool has_key(const key_type& key) const {
    size_t id;
    return get_key_id(key, &id);
  }

  // Get the id of the key in the hash. Returns true if the key is
  // found in the hash, false otherwise.
  inline bool get_key_id(const key_type& key, size_t* id) const {
    key_type        tmp_key;
    const word*     w;
    const offset_t* o;
    return get_key_id(key, id, tmp_key, &w, &o);
  }

  // Optimization version where a tmp_key buffer is provided instead
  // of being allocated. May be faster if many calls to get_key_id are
  // made consecutively by passing the same tmp_key each time.
  inline bool get_key_id(const key_type& key, size_t* id, key_type& tmp_key) const {
    const word*     w;
    const offset_t* o;
    return get_key_id(key, id, tmp_key, &w, &o);
  }

protected:
  // Information and methods to manage the prefetched data.
  struct prefetch_info {
    size_t          id;
    const word*     w;
    const offset_t *o, *lo;
  };
  typedef simple_circular_buffer::pre_alloc<prefetch_info, 8> prefetch_buffer;

  void warm_up_cache(prefetch_buffer& buffer, size_t oid) const {
    buffer.clear();
    for(int i = 0; i < buffer.capacity(); ++i) {
      buffer.push_back();
      prefetch_info& info = buffer.back();
      info.id             = (oid + (i > 0 ? reprobes_[i] : 0)) & size_mask_;
      info.w              = offsets_.word_offset(info.id, &info.o, &info.lo, data_);
      __builtin_prefetch(info.w + info.o->key.woff, 0, 1);
      __builtin_prefetch(info.o, 0, 3);
    }
  }

  void prefetch_next(prefetch_buffer& buffer, size_t oid, uint_t reprobe) const {
    buffer.pop_front();
    //    if(reprobe + buffer.capacity() <= reprobe_limit_.val()) {
      buffer.push_back();
      prefetch_info& info = buffer.back();
      info.id             = (oid + reprobes_[reprobe + buffer.capacity() - 1]) & size_mask_;
      info.w              = offsets_.word_offset(info.id, &info.o, &info.lo, data_);
      __builtin_prefetch(info.w + info.o->key.woff, 0, 1);
      __builtin_prefetch(info.o, 0, 3);
      //    }
  }

public:
  // Optimization version again. Also return the word and the offset
  // information where the key was found. These can be used later one
  // to fetch the value associated with the key.
  inline bool get_key_id(const key_type& key, size_t* id, key_type& tmp_key, const word** w, const offset_t** o) const {
    return get_key_id(key, id, tmp_key, w, o, hash_matrix_.times(key) & size_mask_);
  }

  // Find the actual id of the key in the hash, starting at oid.
  bool get_key_id(const key_type& key, size_t* id, key_type& tmp_key, const word** w, const offset_t** o, const size_t oid) const {
    // This static_assert makes clang++ happy
    static_assert(std::is_pod<prefetch_info>::value, "prefetch_info must be a POD");
    prefetch_info info_ary[prefetch_buffer::capacity()];
    prefetch_buffer buffer(info_ary);
    warm_up_cache(buffer, oid);

    for(uint_t reprobe = 0; reprobe <= reprobe_limit_.val(); ++reprobe) {
      prefetch_info& info = buffer.front();
      key_status st       = get_key_at_id(info.id, tmp_key, info.w, info.o);

      switch(st) {
      case EMPTY:
        return false;
      case FILLED:
        if(oid != tmp_key.get_bits(0, lsize_))
          break;
        tmp_key.template set_bits<false>(0, lsize_, key.get_bits(0, lsize_));
        if(tmp_key != key)
          break;
        *id = info.id;
        *w  = info.w;
        *o  = info.o;
        return true;
      default:
        break;
      }

      prefetch_next(buffer, oid, reprobe + 1);
    } // for

    return false;
  }

  //////////////////////////////
  // Iterator
  //////////////////////////////
  const_iterator begin() { return const_iterator(this); }
  const_iterator begin() const { return const_iterator(this); }
  const_iterator end() { return const_iterator(); }
  const_iterator end() const { return const_iterator(); }

/// Get a slice of an array as an iterator
  template<typename Iterator>
  Iterator iterator_slice(size_t index, size_t nb_slices) const {
    std::pair<size_t, size_t> res = slice(index, nb_slices, size());
    return Iterator(this, res.first, res.second);
  }

  template<typename Iterator>
  Iterator iterator_all() const { return iterator_slice<Iterator>(0, 1); }

  // See hash_counter.hpp for why we added this method. It should not
  // be needed, but I can't get the thing to compile without :(.
  eager_iterator eager_slice(size_t index, size_t nb_slices) const {
    return iterator_slice<eager_iterator>(index, nb_slices);
  }
  region_iterator region_slice(size_t index, size_t nb_slices) const {
    return iterator_slice<region_iterator>(index, nb_slices);
  }

  // Claim a key with the large bit not set. I.e. first entry for a key.
  //
  // id is input/output. Equal to hash & size_maks on input. Equal to
  // actual id where key was set on output. key is already hash
  // shifted and masked to get higher bits. (>> lsize & key_mask)
  // is_new is set on output to true if key did not exists in hash
  // before. *ao points to the actual offsets object and w to the word
  // holding the value.
  bool claim_key(const key_type& key, bool* is_new, size_t* id, const offset_t** _ao, word** _w) {
    uint_t	    reprobe        = 0;
    const offset_t *o, *lo;
    word	   *w, *kw, nkey;
    bool	    key_claimed    = false;
    size_t	    cid            = *id;

    // Akey contains first word of what to store in the key
    // field. I.e. part of the original key (the rest is encoded in
    // the original position) and the reprobe value to substract from
    // the actual position to get to the original position.
    //
    //    MSB                     LSB
    //   +--------------+-------------+
    //   |  MSB of key  |  reprobe    |
    //   + -------------+-------------+
    //     raw_key_len    reprobe_len
    //
    // Akey is updated at every operation to reflect the current
    // reprobe value. nkey is the temporary word containing the part
    // to be stored in the current word kw (+ some offset).
    word      akey          = 1; // start reprobe value == 0. Store reprobe value + 1
    const int to_copy       = std::min((uint16_t)(wsize - offsets_.reprobe_len()), raw_key_len_);
    const int implied_copy  = std::min(key_len_, lsize_);
    akey                   |= key.get_bits(implied_copy, to_copy) << offsets_.reprobe_len();
    const int abits_copied  = implied_copy + to_copy; // Bits from original key already copied, explicitly or implicitly

    do {
      int bits_copied = abits_copied;

      w  = offsets_.word_offset(cid, &o, &lo, data_);
      kw = w + o->key.woff;

      if(o->key.sb_mask1) { // key split on multiple words
        nkey = akey << o->key.boff;
        nkey |= o->key.sb_mask1;
        nkey &= o->key.mask1;

        key_claimed = set_key(kw, nkey, o->key.mask1, o->key.mask1, is_new);
        if(key_claimed) {
          nkey = akey >> o->key.shift;
          if(o->key.full_words) {
            // Copy full words. First one is special
            nkey                  |= key.get_bits(bits_copied, o->key.shift - 1) << (wsize - o->key.shift);
            bits_copied           += o->key.shift - 1;
            nkey                  |= o->key.sb_mask1; // Set bit is MSB
            int copied_full_words  = 1;
            key_claimed            = set_key(kw + copied_full_words, nkey, fmask, fmask, is_new);
            // Copy more full words if needed
            while(bits_copied + wsize - 1 <= key_len_ && key_claimed) {
              nkey               = key.get_bits(bits_copied, wsize - 1);
              bits_copied       += wsize - 1;
              nkey              |= o->key.sb_mask1;
              copied_full_words += 1;
              key_claimed        = set_key(kw + copied_full_words, nkey, fmask, fmask, is_new);
            }
            assert(!key_claimed || (bits_copied < key_len_) == (o->key.sb_mask2 != 0));
            if(o->key.sb_mask2 && key_claimed) { // Copy last word
              nkey               = key.get_bits(bits_copied, key_len_ - bits_copied);
              nkey              |= o->key.sb_mask2;
              copied_full_words += 1;
              key_claimed        = set_key(kw + copied_full_words, nkey, o->key.mask2, o->key.mask2, is_new);
            }
          } else if(o->key.sb_mask2) { // if bits_copied + wsize - 1 < key_len
            // Copy last word, no full words copied
            nkey        |= key.get_bits(bits_copied, key_len_ - bits_copied) << (wsize - o->key.shift);
            nkey        |= o->key.sb_mask2;
            nkey        &= o->key.mask2;
            key_claimed  = set_key(kw + 1, nkey, o->key.mask2, o->key.mask2, is_new);
          }
        } // if(key_claimed)
      } else { // key on one word
        nkey = akey << o->key.boff;
        nkey &= o->key.mask1;
        key_claimed = set_key(kw, nkey, o->key.mask1, o->key.mask1, is_new);
      }
      if(!key_claimed) { // reprobe
        if(++reprobe > reprobe_limit_.val())
          return false;
        cid = (*id + reprobes_[reprobe]) & size_mask_;
        akey = (akey & ~offsets_.reprobe_mask()) | (reprobe + 1);
      }
    } while(!key_claimed);

    *id  = cid;
    *_w  = w;
    *_ao = o;
    return true;
  }

  // Claim large key. Enter an entry for a key when it is not the
  // first entry. Only encode the number of reprobe hops back to the
  // first entry of the key in the hash table. It is simpler as can
  // takes less than one word in length.
  bool claim_large_key(size_t* id, const offset_t** _ao, word** _w) {
    uint_t          reprobe     = 0;
    size_t          cid         = *id;
    const offset_t *o, *lo;
    word           *w, *kw, nkey;
    bool            key_claimed = false;

    do {
      w = offsets_.word_offset(cid, &o, &lo, data_);
      kw = w + lo->key.woff;

      if(lo->key.sb_mask1) { // key split on multiple words
        nkey = (reprobe << lo->key.boff) | lo->key.sb_mask1 | lo->key.lb_mask;
        nkey &= lo->key.mask1;

        // Use o->key.mask1 and not lo->key.mask1 as the first one is
        // guaranteed to be bigger. The key needs to be free on its
        // longer mask to claim it!
        key_claimed = set_key(kw, nkey, o->key.mask1, lo->key.mask1);
        if(key_claimed) {
          nkey         = (reprobe >> lo->key.shift) | lo->key.sb_mask2;
          nkey        &= lo->key.mask2;
          key_claimed  = set_key(kw + 1, nkey, o->key.full_words ? fmask : o->key.mask2, lo->key.mask2);
        }
      } else { // key on 1 word
        nkey  = (reprobe << lo->key.boff) | lo->key.lb_mask;
        nkey &= lo->key.mask1;
        key_claimed = set_key(kw, nkey, o->key.mask1, lo->key.mask1);
      }
      if(!key_claimed) { //reprobe
        if(++reprobe > reprobe_limit_.val())
          return false;
        cid  = (*id + reprobes_[reprobe]) & size_mask_;
      }
    } while(!key_claimed);

    *id  = cid;
    *_w  = w;
    *_ao = lo;
    return true;
  }

  // Add val to key. id is the starting place (result of hash
  // computation). eid is set to the effective place in the
  // array. large is set to true is setting a large key (upon
  // recurrence if there is a carry).
  bool add_rec(size_t id, const key_type& key, word val, bool large, bool* is_new, size_t* eid, unsigned int* carry_shift) {
    const offset_t *ao = 0;
    word	   *w  = 0;

    bool claimed = false;
    if(large)
      claimed = claim_large_key(&id, &ao, &w);
    else
      claimed = claim_key(key, is_new, &id, &ao, &w);
    if(!claimed)
      return false;
    *eid = id;
    return add_rec_at(id, key, val, ao, w, carry_shift);
  }

  bool add_rec_at(size_t id, const key_type& key, word val, const offset_t* ao, word* w, unsigned int* carry_shift) {
    // Increment value
    word *vw       = w + ao->val.woff;
    word  cary     = add_val(vw, val, ao->val.boff, ao->val.mask1);
    cary         >>= ao->val.shift;
    *carry_shift  += ao->val.shift;
    if(cary && ao->val.mask2) { // value split on two words
      cary           = add_val(vw + 1, cary, 0, ao->val.mask2);
      cary         >>= ao->val.cshift;
      *carry_shift  += ao->val.cshift;
    }
    if(!cary)
      return true;

    id = (id + reprobes_[0]) & size_mask_;
    size_t ignore_eid;
    bool   ignore_is_new;
    return add_rec(id, key, cary, true, &ignore_is_new, &ignore_eid, carry_shift);

      // // Adding failed, table is full. Need to back-track and
      // // substract val.
      //      std::cerr << "Failed to add large part of value -> return false\n";
      // cary = add_val(vw, ((word)1 << offsets_.val_len()) - val,
      //                ao->val.boff, ao->val.mask1);
      // cary >>= ao->val.shift;
      // if(cary && ao->val.mask2) {
      //   // Can I ignore the cary here? Table is known to be full, so
      //   // not much of a choice. But does it leave the table in a
      //   // consistent state?
      //   add_val(vw + 1, cary, 0, ao->val.mask2);
      // }
      //      return false;
  }

  // Atomic methods to set the key. Attempt to set nkey in word w. All
  // bits matching free_mask must be unset and the bits matching
  // equal_mask must be equal for a success in setting the key. Set
  // is_new to true if the spot was previously empty. Otherwise, if
  // is_new is false but true is returned, the key was already present
  // at that spot.
  inline bool set_key(word *w, word nkey, word free_mask, word equal_mask, bool *is_new) {
    word ow = *w, nw, okey;

    okey = ow & free_mask;
    while(okey == 0) { // large bit not set && key is free
      nw = atomic_.cas(w, ow, ow | nkey);
      if(nw == ow) {
        *is_new = true;
        return true;
      }
      ow = nw;
      okey = ow & free_mask;
    }
    *is_new = false;
    return (ow & equal_mask) == nkey;
  }

  inline bool set_key(word *w, word nkey, word free_mask, word equal_mask) {
    bool is_new;
    return set_key(w, nkey, free_mask, equal_mask, &is_new);
  }

  // Add val the value in word w, with shift and mask giving the
  // particular part of the word in which the value is stored. The
  // return value is the carry.
  inline word add_val(word *w, word val, uint_t shift, word mask) {
    word now = *w, ow, nw, nval;

    do {
      ow = now;
      nval = ((ow & mask) >> shift) + val;
      nw = (ow & ~mask) | ((nval << shift) & mask);
      now = atomic_.cas(w, ow, nw);
    } while(now != ow);

    return nval & (~(mask >> shift));
  }

  // Return the key and value at position id. If the slot at id is
  // empty or has the large bit set, returns false. Otherwise, returns
  // the key and the value is the sum of all the entries in the hash
  // table for that key. I.e., the table is search forward for entries
  // with large bit set pointing back to the key at id, and all those
  // values are summed up.
  key_status get_key_val_at_id(size_t id, key_type& key, word& val, const bool carry_bit = false) const {
    const word*     w;
    const offset_t* o;

    key_status st = get_key_at_id(id, key, &w, &o);
    if(st != FILLED)
       return st;

    val = get_val_at_id(id, w, o, true, carry_bit);

    return FILLED;
  }

  // Get a the key at the given id. It also returns the word and
  // offset information in w and o. The return value is EMPTY (no key
  // at id), FILLED (there is a key at id), LBSET (the large bit is
  // set, hence the key is only a pointer back to the real key).
  //
  // The key returned contains the original id in the hash as its
  // lsize_ lsb bits. To obtain the full key, one needs to compute the
  // product with the inverse matrix to get the lsb bits.
  inline key_status get_key_at_id(size_t id, key_type& key, const word** w, const offset_t** o) const {
    const offset_t *lo;
    *w = offsets_.word_offset(id, o, &lo, data_);
    return get_key_at_id(id, key, *w, *o);
  }

  // Sam as above, but it assume that the word w and o for id have
  // already be computed (like already prefetched).
  key_status get_key_at_id(size_t id, key_type&key, const word* w, const offset_t* o) const {
    const word*     kvw      = w + o->key.woff;
    word            key_word = *kvw;
    word            kreprobe = 0;

    const key_offsets& key_o = o->key;
    if(key_word & key_o.lb_mask)
      return LBSET;
    const int implied_copy = std::min(lsize_, key_len_);
    int       bits_copied  = implied_copy;
    if(key_o.sb_mask1) {
      if((key_word & key_o.sb_mask1) == 0)
        return EMPTY;
      kreprobe = (key_word & key_o.mask1 & ~key_o.sb_mask1) >> key_o.boff;
      if(key_o.full_words) {
        // Copy full words. First one is special
        key_word = *(kvw + 1);
        if(offsets_.reprobe_len() < key_o.shift) {
          key.set_bits(bits_copied, key_o.shift - offsets_.reprobe_len(), kreprobe >> offsets_.reprobe_len());
          bits_copied += key_o.shift - offsets_.reprobe_len();
          kreprobe    &= offsets_.reprobe_mask();
          key.set_bits(bits_copied, wsize - 1, key_word & ~key_o.sb_mask1);
          bits_copied += wsize - 1;
        } else {
          int reprobe_left  = offsets_.reprobe_len() - key_o.shift;
          kreprobe         |= (key_word & (((word)1 << reprobe_left) - 1)) << key_o.shift;
          key.set_bits(bits_copied, wsize - 1 - reprobe_left, (key_word & ~key_o.sb_mask1) >> reprobe_left);
          bits_copied += wsize - 1 - reprobe_left;
        }
        int word_copied = 2;
        while(bits_copied + wsize - 1 <= key_len_) {
          key.set_bits(bits_copied, wsize - 1, *(kvw + word_copied++) & (fmask >> 1));
          bits_copied += wsize - 1;
        }
        if(key_o.sb_mask2)
          key.set_bits(bits_copied, key_len_ - bits_copied, *(kvw + word_copied) & key_o.mask2 & ~key_o.sb_mask2);
      } else if(key_o.sb_mask2) { // if(bits_copied + wsize - 1 < key_len
        // Two words but no full words
        key_word = *(kvw + 1) & key_o.mask2 & ~key_o.sb_mask2;
        if(offsets_.reprobe_len() < key_o.shift) {
          key.set_bits(bits_copied, key_o.shift - offsets_.reprobe_len(), kreprobe >> offsets_.reprobe_len());
          bits_copied += key_o.shift - offsets_.reprobe_len();
          kreprobe    &= offsets_.reprobe_mask();
          key.set_bits(bits_copied, key_len_ - bits_copied, key_word);
        } else {
          int reprobe_left  = offsets_.reprobe_len() - key_o.shift;
          kreprobe         |= (key_word & (((word)1 << reprobe_left) - 1)) << key_o.shift;
          key.set_bits(bits_copied, key_len_ - bits_copied, key_word >> reprobe_left);
        }
      }
    } else { // if(key_o.sb_mask1
      // Everything in 1 word
      key_word = (key_word & key_o.mask1) >> key_o.boff;
      if(key_word == 0)
        return EMPTY;
      kreprobe = key_word & offsets_.reprobe_mask();
      key.set_bits(bits_copied, raw_key_len_, key_word >> offsets_.reprobe_len());
    }
    // Compute missing oid so that the original key can be computed
    // back through the inverse matrix. Although the key may have a
    // length of key_len_, which may be less than lsize_, assume that
    // it still fit here as lsize_ is less than a word length. Need all lsize_.
    size_t oid = id; // Original id
    if(kreprobe > 1)
      oid -= reprobes_[kreprobe - 1];
    oid &= size_mask_;
    // Can use more bits than mer size. That's OK, will fix it later
    // when computing the actual mers by computing the product with
    // the inverse matrix.
    key.template set_bits<0>(0, lsize_, oid);

    return FILLED;
  }

  word get_val_at_id(const size_t id, const word* w, const offset_t* o, const bool reprobe = true,
                     const bool carry_bit = false) const {
    word            val = 0;
    if(val_len() == 0)
      return val;

    // First part of value
    const word* kvw = w + o->val.woff;
    val = ((*kvw) & o->val.mask1) >> o->val.boff;
    if(o->val.mask2)
      val |= ((*(kvw+1)) & o->val.mask2) << o->val.shift;

    // Do we want to get the large value
    bool do_reprobe = reprobe;
    if(carry_bit && do_reprobe) {
      do_reprobe   = do_reprobe && (val & 0x1);
      val        >>= 1;
    }
    if(!do_reprobe)
      return val;

    return resolve_val_rec((id + reprobes_[0]) & size_mask_, val, carry_bit);
  }

  word resolve_val_rec(const size_t id, word val, const bool carry_bit, const uint_t overflows = 0) const {
    uint_t          reprobe = 0;
    size_t          cid     = id;

    while(reprobe <= reprobe_limit_.val()) {
      const offset_t     *o, *lo;
      const word*         w    = offsets_.word_offset(cid, &o, &lo, data_);
      const word*         kw   = w + o->key.woff;
      word                nkey = *kw;
      const key_offsets&  lkey = lo->key;

      if(nkey & lkey.lb_mask) {
        // If the large bit is set, the size of the key (reprobe_len)
        // is guaranteed to have a length of at most 1 word.
        if(lkey.sb_mask1) {
          nkey  = (nkey & lkey.mask1 & ~lkey.sb_mask1) >> lkey.boff;
          nkey |= ((*(kw+1)) & lkey.mask2 & ~lkey.sb_mask2) << lkey.shift;
        } else {
          nkey = (nkey & lkey.mask1) >> lkey.boff;
        }
        if(nkey == reprobe) {
          const val_offsets& lval = lo->val;
          const word*        vw   = w + lval.woff;
          word               nval = ((*vw) & lval.mask1) >> lval.boff;
          if(lval.mask2)
            nval |= ((*(vw+1)) & lval.mask2) << lval.shift;

          bool do_reprobe = true;
          if(carry_bit) {
            do_reprobe   = nval & 0x1;
            nval       >>= 1;
          }

          nval <<= offsets_.val_len();
          nval <<= offsets_.lval_len() * overflows;
          val   += nval;

          if(!do_reprobe)
            return val;

          return resolve_val_rec((cid + reprobes_[0]) & size_mask_, val, carry_bit, overflows + 1);
        }
      } else if((nkey & o->key.mask1) == 0) {
        break;
      }

      cid  = (id + reprobes_[++reprobe]) & size_mask_;
    }

    return val;
  }

};

template<typename Key, typename word = uint64_t, typename atomic_t = ::atomic::gcc, typename mem_block_t = ::allocators::mmap>
class array :
    protected mem_block_t,
    public array_base<Key, word, atomic_t, array<Key, word, atomic_t, mem_block_t> >
{
  typedef array_base<Key, word, atomic_t, array<Key, word, atomic_t, mem_block_t> > super;
  friend class array_base<Key, word, atomic_t, array<Key, word, atomic_t, mem_block_t> >;

public:
  array(size_t size, // Size of hash. To be rounded up to a power of 2
        uint16_t key_len, // Size of key in bits
        uint16_t val_len, // Size of val in bits
        uint16_t reprobe_limit, // Maximum reprobe
        const size_t* reprobes = quadratic_reprobes) : // Reprobing policy
    mem_block_t(),
    super(size, key_len, val_len, reprobe_limit, RectangularBinaryMatrix(ceilLog2(size), key_len).randomize_pseudo_inverse(),
          reprobes)
  { }

protected:
  word* alloc_data(size_t s) {
    mem_block_t::realloc(s);
    return (word*)mem_block_t::get_ptr();
  }
};

struct ptr_info {
  void*  ptr_;
  size_t bytes_;
  ptr_info(void* ptr, size_t bytes) : ptr_(ptr), bytes_(bytes) { }
};
template<typename Key, typename word = uint64_t, typename atomic_t = ::atomic::gcc>
class array_raw :
    protected ptr_info,
    public array_base<Key, word, atomic_t, array<Key, word, atomic_t> >
{
  typedef array_base<Key, word, atomic_t, array<Key, word, atomic_t> > super;
  friend class array_base<Key, word, atomic_t, array<Key, word, atomic_t> >;

public:
  array_raw(void* ptr,
            size_t bytes, // Memory available at ptr
            size_t size, // Size of hash in number of entries. To be rounded up to a power of 2
            uint16_t key_len, // Size of key in bits
            uint16_t val_len, // Size of val in bits
            uint16_t reprobe_limit, // Maximum reprobe
            RectangularBinaryMatrix m,
            const size_t* reprobes = quadratic_reprobes) : // Reprobing policy
    ptr_info(ptr, bytes),
    super(size, key_len, val_len, reprobe_limit, m, reprobes)
  { }

protected:
  word* alloc_data(size_t s) {
    assert(bytes_ == s);
    return (word*)ptr_;
  }
};

} } // namespace jellyfish { namespace large_hash_array

#endif /* __JELLYFISH_LARGE_HASH_ARRAY_HPP__ */
