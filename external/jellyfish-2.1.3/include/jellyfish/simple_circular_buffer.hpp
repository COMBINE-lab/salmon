/* Jellyfish
 * Copyright (C) 2012  Genome group at University of Maryland.
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __SIMPLE_CIRCULAR_BUFFER_H__
#define __SIMPLE_CIRCULAR_BUFFER_H__

#include <memory>

namespace jellyfish {
  namespace simple_circular_buffer {

    // T: type of element in container. D: type of derived class for
    // CRTP. A: allocator type.
    template<typename T, typename D>
    class base {
    public:
      explicit base(T* data) :
        data_(data), front_(0), back_(0), full_(false)
      { }

      // Return true if empty
      bool empty() const {
        return front_ == back_ && !full();
      }
      // Return true if full
      bool full() const {
        return full_;
      }
      void clear() {
        front_ = back_;
        full_  = false;
      }

      // Valid only if empty() is false
      T& front() {
        return data_[front_];
      }
      // Valid only if empty() is false
      T& back() {
        return data_[prev_index(back_)];
      }

      // Unlike the corresponding method on list or deqeue, push_back may
      // fail if full() is true. Then false is returned.
      bool push_back(const T& x) {
        if(full())
          return false;
        data_[back_] = x;
        back_        = next_index(back_);
        full_        = back_ == front_;
        return true;
      }

      bool push_back() {
        if(full())
          return false;
        back_ = next_index(back_);
        full_ = back_ == front_;
        return true;
      }

      // Pop an element from the front. It has no effect if empty() is true
      void pop_front() {
        if(empty())
          return;
        front_ = next_index(front_);
        full_  = false;
      }

      int size() const {
        if(full())
          return static_cast<const D*>(this)->capacity();
        int s = back_ - front_;
        return s < 0 ? s + static_cast<const D*>(this)->capacity() : s;
      }

    protected:
      int next_index(int i) const {
        return (i + 1) % static_cast<const D*>(this)->capacity();
      }
      int prev_index(int i) const {
        return i ? i - 1 : static_cast<const D*>(this)->capacity() - 1;
      }
      T* data() const { return data_; }

      T*   data_;
      int  front_, back_;
      bool full_;
    };

    template<typename T, int capa>
    class pre_alloc : public base<T, pre_alloc<T, capa> > {
      typedef base<T, pre_alloc<T, capa> > super;
    public:
      explicit pre_alloc(T* data) : super(data) { }
      static int capacity() { return capa; }
    };

    // template<typename T, int capa, typename A = std::allocator<T> >
    // class fixed : public base<T, fixed<T, capa, A>, A> {
    //   typedef base<T, fixed<T, capa, A>, A> super;
    // public:
    //   explicit fixed(const T v = T()) : super(capa, v) { }
    //   //    fixed(const int ignored_size, const T v = T()) : super(capa, v) { }

    //   int capacity() const { return capa; }
    // };

    // template<typename T, typename A = std::allocator<T> >
    // class dyn : public base<T, dyn<T, A>, A> {
    //   typedef base<T, dyn<T, A>, A> super;
    // public:
    //   explicit dyn(int size, const T v = T()) : super(size, v), capa_(size) { }

    //   int capacity() const { return capa_; }
    //   int capa_;
    // };
  }
}
#endif /* __SIMPLE_CIRCULAR_BUFFER_H__ */
