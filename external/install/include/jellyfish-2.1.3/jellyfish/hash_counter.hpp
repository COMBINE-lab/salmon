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

#ifndef __HASH_COUNTER_HPP__
#define __HASH_COUNTER_HPP__

#include <stdexcept>

#include <jellyfish/large_hash_array.hpp>
#include <jellyfish/locks_pthread.hpp>
#include <jellyfish/dumper.hpp>

/// Cooperative version of the hash_counter. In this implementation,
/// it is expected that the given number of threads will call the
/// `add` method regularly. In case the hash table is full, it gets
/// enlarged using all the threads. After the work is done, every
/// thread must call promptly the `done` method.

namespace jellyfish{ namespace cooperative {

template<typename Key, typename word = uint64_t, typename atomic_t = ::atomic::gcc, typename mem_block_t = ::allocators::mmap>
class hash_counter {
public:
  typedef typename large_hash::array<Key, word, atomic_t, mem_block_t> array;
  typedef typename array::key_type                                     key_type;
  typedef typename array::mapped_type                                  mapped_type;
  typedef typename array::value_type                                   value_type;
  typedef typename array::reference                                    reference;
  typedef typename array::const_reference                              const_reference;
  typedef typename array::pointer                                      pointer;
  typedef typename array::const_pointer                                const_pointer;
  typedef typename array::eager_iterator                               eager_iterator;
  typedef typename array::lazy_iterator                                lazy_iterator;

protected:
  array*                  ary_;
  array*                  new_ary_;
  uint16_t                nb_threads_;
  locks::pthread::barrier size_barrier_;
  volatile uint16_t       size_thid_, done_threads_;
  bool                    do_size_doubling_;
  dumper_t<array>*        dumper_;

public:
  hash_counter(size_t size, // Size of hash. To be rounded up to a power of 2
               uint16_t key_len, // Size of key in bits
               uint16_t val_len, // Size of val in bits
               uint16_t nb_threads, // Number of threads accessing this hash
               uint16_t reprobe_limit = 126, // Maximum reprobe
               const size_t* reprobes = jellyfish::quadratic_reprobes) :
    ary_(new array(size, key_len, val_len, reprobe_limit, reprobes)),
    new_ary_(0),
    nb_threads_(nb_threads),
    size_barrier_(nb_threads),
    size_thid_(0),
    done_threads_(0),
    do_size_doubling_(true),
    dumper_(0)
  { }

  ~hash_counter() {
    delete ary_;
  }

  array* ary() { return ary_; }
  const array* ary() const { return ary_; }
  size_t size() const { return ary_->size(); }
  uint16_t key_len() const { return ary_->key_len(); }
  uint16_t val_len() const { return ary_->val_len(); }
  uint16_t nb_threads() const { return nb_threads; }
  uint16_t reprobe_limit() const { return ary_->max_reprobe(); }


  /// Whether we attempt to double the size of the hash when full.
  bool do_size_doubling() const { return do_size_doubling_; }
  /// Set whether we attempt to double the size of the hash when full.
  void do_size_doubling(bool v) { do_size_doubling_ = v; }

  /// Set dumper responsible for cleaning out the array.
  void dumper(dumper_t<array> *d) { dumper_ = d; }

  /// Add `v` to the entry `k`. It returns in `is_new` true if the
  /// entry `k` did not exist in the hash. In `id` is returned the
  /// final position of `k` in the hash array.
  void add(const Key& k, uint64_t v, bool* is_new, size_t* id) {
    unsigned int carry_shift  = 0;
    bool*        is_new_ptr   = is_new;
    size_t*      id_ptr       = id;
    bool         is_new_void  = false;
    size_t       id_void      = false;

    while(!ary_->add(k, v, &carry_shift, is_new_ptr, id_ptr)) {
      handle_full_ary();
      v &= ~(uint64_t)0 << carry_shift;
      // If carry_shift == 0, failed to allocate the first field for
      // key, hence status of is_new and value for id are not
      // determined yet. On the other hand, if carry_shift > 0, we
      // failed while adding extra field for large key, so the status
      // of is_new and value of id are known. We do not update them in future
      // calls.
      if(carry_shift) {
        is_new_ptr = &is_new_void;
        id_ptr     = &id_void;
      }
    }
  }

  /// Add `v` to the entry `k`. This method is multi-thread safe. If
  /// the entry for `k` does not exists, it is inserted.
  ///
  /// @param k Key to add to
  /// @param v Value to add
  inline void add(const Key& k, uint64_t v) {
    bool   is_new;
    size_t id;
    add(k, v, &is_new, &id);
  }

  /// Insert the key `k` in the hash. The value is not changed or set
  /// to 0 if not already in the hash.
  ///
  /// @param k Key to insert
  inline void set(const Key& k) {
    bool   is_new;
    size_t id;
    set(k, &is_new, &id);
  }

  /// Insert the key `k` in the hash. The value is not changed or set
  /// to 0 if not already in the hash. Set `is_new` to true if `k` did
  /// not already exist in the hash. In `id` is returned the final
  /// position of `k` in the hash.
  void set(const Key& k, bool* is_new, size_t* id) {
    while(!ary_->set(k, is_new, id))
      handle_full_ary();
  }

  /// Update the value of key `k` by adding `v`, if `k` is already
  /// present in the hash, otherwise this nothing happens. Returns
  /// true if `k` is already in the hash, false otherwise.
  bool update_add(const Key& k, uint64_t v) {
    Key tmp_key;
    return update_add(k, v, tmp_key);
  }

  bool update_add(const Key& k, uint64_t v, Key& tmp_key) {
    unsigned int carry_shift = 0;

    while(true) {
      if(ary_->update_add(k, v, &carry_shift, tmp_key))
        return true;
      if(carry_shift == 0)
        return false;
      handle_full_ary();
      v &= ~(uint64_t)0 << carry_shift;
    }
  }

  /// Signify that thread is done and wait for all threads to be done.
  void done() {
    atomic_t::fetch_add(&done_threads_, (uint16_t)1);
    while(!handle_full_ary()) ;
  }

protected:
  // Double the size of the hash and return false. Unless all the
  // thread have reported they are done, in which case do nothing and
  // return true.
  bool handle_full_ary() {
    bool serial_thread = size_barrier_.wait();
    if(done_threads_ >= nb_threads_) // All done?
      return true;

    bool success = false;
    if(do_size_doubling_)
      success = success || double_size(serial_thread);

    if(!success && dumper_) {
      if(serial_thread)
        dumper_->dump(ary_);
      success = true;
      size_barrier_.wait();
    }

    if(!success)
      throw std::runtime_error("Hash full");

    return false;
  }

  bool double_size(bool serial_thread) {
    if(serial_thread) {// Allocate new array for size doubling
      try {
        new_ary_   = new array(ary_->size() * 2, ary_->key_len(), ary_->val_len(),
                               ary_->max_reprobe(), ary_->reprobes());
       } catch(typename array::ErrorAllocation e) {
        new_ary_ = 0;
      }
    }
    size_thid_ = 0;

    size_barrier_.wait();
    array* my_ary = *(array* volatile*)&new_ary_;
    if(!my_ary) // Allocation failed
      return false;

    // Copy data from old to new
    uint16_t       id = atomic_t::fetch_add(&size_thid_, (uint16_t)1);
    // Why doesn't the following work? Seems like a bug to
    // me. Equivalent call works in test_large_hash_array. Or am I
    // missing something?
    // eager_iterator it = ary_->iterator_slice<eager_iterator>(id, nb_threads_);
    eager_iterator it = ary_->eager_slice(id, nb_threads_);
    while(it.next())
      my_ary->add(it.key(), it.val());

    size_barrier_.wait();

    if(serial_thread) { // Set new ary to be current and free old
      delete ary_;
      ary_ = new_ary_;
    }

    // Done. Last sync point
    size_barrier_.wait();
    return true;
  }
};

} } // namespace jellyfish { namespace cooperative {
#endif /* __HASH_COUNTER_HPP__ */
