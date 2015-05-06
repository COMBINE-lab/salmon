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

#ifndef __JELLYFISH_CPP_ARRAY_HPP_
#define __JELLYFISH_CPP_ARRAY_HPP_

#include <cstddef>
#include <memory>

namespace jellyfish {

/// Fix length array of type T. An element is initialized with the init method.
///   new (this->data() + i) T(
template<typename T>
class cpp_array {
protected:
  std::pair<T*, ptrdiff_t> data_;
  std::pair<bool*, ptrdiff_t> init_;
  size_t size_;

public:
  cpp_array(size_t size) :
  data_(std::get_temporary_buffer<T>(size)),
  init_(std::get_temporary_buffer<bool>(size)),
  size_(size) {
    if(data_.first == 0 || init_.first == 0) {
      std::return_temporary_buffer(data_.first);
      std::return_temporary_buffer(init_.first);
      throw std::bad_alloc();
    }
    memset(init_.first, '\0', sizeof(bool) * size_);
  }

  ~cpp_array() {
    clear();
    std::return_temporary_buffer(data_.first);
    std::return_temporary_buffer(init_.first);
  }

  /// Initialize element i with 0 argument
  void init(size_t i) {
    release(i);
    new (data_.first + i) T();
    init_.first[i] = true;
  }

  /// Initialize element i with 1 argument
  template<typename A1>
  void init(size_t i, A1& a1) {
    release(i);
    new (data_.first + i) T(a1);
    init_.first[i] = true;
  }
  template<typename A1>
  void init(size_t i, A1* a1) {
    release(i);
    new (data_.first + i) T(a1);
    init_.first[i] = true;
  }
  /// Initialize element i with 2 arguments
  template<typename A1, typename A2>
  void init(size_t i, A1& a1, A2& a2) {
    release(i);
    new (data_.first + i) T(a1, a2);
    init_.first[i] = true;
  }
  template<typename A1, typename A2>
  void init(size_t i, A1* a1, A2& a2) {
    release(i);
    new (data_.first + i) T(a1, a2);
    init_.first[i] = true;
  }
  template<typename A1, typename A2>
  void init(size_t i, A1& a1, A2* a2) {
    release(i);
    new (data_.first + i) T(a1, a2);
    init_.first[i] = true;
  }
  template<typename A1, typename A2>
  void init(size_t i, A1* a1, A2* a2) {
    release(i);
    new (data_.first + i) T(a1, a2);
    init_.first[i] = true;
  }

  /// Initialize element i with 3 arguments
  template<typename A1, typename A2, typename A3>
  void init(size_t i, A1 a1, A2 a2, A3 a3) {
    release(i);
    new (data_.first + i) T(a1, a2, a3);
    init_.first[i] = true;
  }
  /// Initialize element i with 4 arguments
  template<typename A1, typename A2, typename A3, typename A4>
  void init(size_t i, A1 a1, A2 a2, A3 a3, A4 a4) {
    release(i);
    new (data_.first + i) T(a1, a2, a3, a4);
    init_.first[i] = true;
  }
  /// Initialize element i with 5 arguments
  template<typename A1, typename A2, typename A3, typename A4, typename A5>
  void init(size_t i, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5) {
    release(i);
    new (data_.first + i) T(a1, a2, a3, a4, a5);
    init_.first[i] = true;
  }

  void release(size_t i) {
    if(init_.first[i]) {
      data_.first[i].~T();
      init_.first[i] = false;
    }
  }

  size_t size() const { return size_; }
  bool empty() const { return size_ == 0; }
  T& operator[](size_t i) { return data_.first[i]; }
  const T& operator[](size_t i) const { return data_.first[i]; }
  bool initialized(size_t i) const { return init_.first[i]; }

  T* begin() { return data_.first; }
  T* end() { return data_.first + size_; }
  const T* begin() const { return data_.first; }
  const T* end() const { return data_.end + size_; }
  const T* cbegin() const { return data_.first; }
  const T* cend() const { return data_.end + size_; }

  T* data() { return data_.first; }
  const T* data() const { return data_.first; }

  T& front() { return data_.first[0]; }
  T& back() { return data_.first[size_ - 1]; }
  const T& front() const { return data_.first[0]; }
  const T& back() const { return data_.first[size_ - 1]; }

  void clear() {
    for(size_t i = 0; i < size_; ++i)
      release(i);
  }
};
} // namespace jellyfish

#endif /* __JELLYFISH_CPP_ARRAY_HPP_ */
