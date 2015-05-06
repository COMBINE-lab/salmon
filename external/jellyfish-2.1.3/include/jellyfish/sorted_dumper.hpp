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

#ifndef __JELLYFISH_SORTED_DUMPER_HPP__
#define __JELLYFISH_SORTED_DUMPER_HPP__

#include <iostream>
#include <sstream>
#include <fstream>

#include <jellyfish/dumper.hpp>
#include <jellyfish/mer_heap.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/token_ring.hpp>
#include <jellyfish/locks_pthread.hpp>
#include <jellyfish/file_header.hpp>

namespace jellyfish {
/// Sorted dumper. Write mers according to the hash order. It
/// implements the CRTP to effectively write the k-mer/value pair.
template<typename D, typename storage_t>
class sorted_dumper : public dumper_t<storage_t>, public thread_exec {
protected:
  typedef typename storage_t::key_type key_type;
  typedef typename storage_t::region_iterator iterator;
  typedef typename mer_heap::heap<typename    storage_t::key_type, iterator> heap_type;
  typedef typename heap_type::const_item_t    heap_item;
  typedef jellyfish::token_ring<locks::pthread::cond> token_ring;
  typedef typename token_ring::token token_type;

  int                       nb_threads_;
  token_ring                ring_;
  const char*               file_prefix_;
  storage_t*                ary_;
  file_header*              header_;
  bool                      zero_array_;
  std::ofstream             out_;
  std::pair<size_t, size_t> block_info; // { nb blocks, nb records }

public:
  sorted_dumper(int nb_threads, const char* file_prefix, file_header* header = 0) :
    nb_threads_(nb_threads),
    ring_(nb_threads),
    file_prefix_(file_prefix),
    header_(header),
    zero_array_(true)
  { }

  bool zero_array() const { return zero_array_; }
  void zero_array(bool v) { zero_array_ = v; }

  virtual void _dump(storage_t* ary) {
    ary_ = ary;
    block_info = ary_->blocks_for_records(5 * ary_->max_reprobe_offset());

    this->open_next_file(file_prefix_, out_);
    if(header_)
      header_->write(out_);

    ring_.reset();
    exec_join(nb_threads_);
    out_.close();
    if(zero_array_)
      ary_->zero_blocks(0, block_info.first); // zero out last group of blocks
  }

  virtual void start(const int i) {
    std::ostringstream           buffer;
    heap_type                    heap(ary_->max_reprobe_offset());
    token_type&                  token = ring_[i];
    size_t                       count = 0;
    typename storage_t::key_type key;

    for(size_t id = i; id * block_info.second < ary_->size(); id += nb_threads_) {
      // Fill buffer
      iterator it(ary_, id * block_info.second, (id + 1) * block_info.second, key);
      heap.fill(it);

      while(heap.is_not_empty()) {
        heap_item item = heap.head();
        if(item->val_ >= this->min() && item->val_ <= this->max())
          static_cast<D*>(this)->write_key_value_pair(buffer, item);
        ++count;
        heap.pop();
        if(it.next())
          heap.push(it);
      }

      // Write buffer
      token.wait();
      out_.write(buffer.str().data(), buffer.tellp());
      token.pass();

      buffer.seekp(0);
      if(id > 0 && zero_array_)
        ary_->zero_blocks(id * block_info.first, block_info.first);
    }
  }
};
} // namespace jellyfish {

#endif /* __JELLYFISH_SORTED_DUMPER_HPP__ */
