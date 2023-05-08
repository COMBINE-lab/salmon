// itlib-dynamic-bitset v1.02
//
// A class similar to std::bitset but the size is not a part of the type
//
// SPDX-License-Identifier: MIT
// MIT License:
// Copyright(c) 2020-2021 Borislav Stanimirov
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files(the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and / or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions :
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT.IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
//
//                  VERSION HISTORY
//
//  1.02 (2021-01-29) Added copy
//                    Added iterator +-num arithmetic
//  1.01 (2021-01-28) Added assign
//                    Added append
//                    Clear trailing bits after resize
//  1.00 (2020-10-18) Initial release
//
//
//                  DOCUMENTATION
//
// Simply include this file wherever you need.
// It defines the class itlib::dynamic_bitset, which similar to std::bitset,
// but does not have the size in its type. It's also similar to
// std::vector<bool> but has more limited ways for modification.
//
// To construct itlib::dynamic_bitset you need to provide a std::vector-like
// container of unsigned integers as a template parameter. It needs to have
// random access. The value_type of the container determines the word type of
// the dynamic_bitset type.
//
// Methods:
// * dynamic_bitset(const allocator_type& a = {})
//      construct and empty bitset
// * dynamic_bitset(size_type size, word_type value = 0, const allocator_type& alloc = {})
//      construct a bitset of size s and intialize per word
//
// * copy, and move ctors and assignment operators
//
// * size_type size() const noexcept
//      number of bits in buffer
// * buffer_size_type word_size() const noexcept
//      size of underlying buffer
// * word_type* data() noexcept
// * const word_type* data() const noexcept
//      data of buffer
// * bool empty() const noexcept
//      empty check
// * size_t byte_size() const noexcept
//      number of bytes to store buffer. Note that this is not necessarily a
//      multiple of sizeof(word_type)
// * const buffer_type& buffer() const
//      get underlying buffer
// * bool test(size_type i) const noexcept
//      test single bit
// * bool all() const noexcept
//      all bits true?
// * bool any() const noexcept
//      any bits true?
// * bool none() const noexcept
//      no bits true?
// * void set(size_type i)
//      set bit to true
// * void reset(size_type i)
//      set bit to false
// * void set(size_type i, bool b)
//      set bit to value
// * void flip(size_type i)
//      flip bit
// * void assign(const Buffer&)
//      set values from a given word-type buffer
// * void append(const dynamic_bitset&)
//      append an existing dynamic bitset to this one
// * void reserve(size_type size)
//      reserve buffer for size bits
// * void resize(size_type size)
//      resize to size bits
// * void push_back(bool b)
//      add bit and resize
// * bool operator[](size_type i) const
//      get bit at index
// * bitref operator[](size_type i)
//      get bit reference. `bitvec[231] = true;` is allowed
// * iterator begin()
// * iterator end()
// * const_iterator begin() const
// * const_iterator end() const
//      const and mutating iterators
// * void copy(iterator from, iterator to, const_iterator source)
//      copy elements from another dynamic_bitset
//
// * static constexpr buffer_size_type word_size(size_type size) noexcept
//      required word size for a given bit size
// * static constexpr buffer_size_type word_index(size_type i) noexcept
//      index of a word in buffer for a given bit index
// * static constexpr uint8_t bit_index(size_type i) noexcept
//      index of a bit within its word
// * static constexpr word_type word_mask(size_type index) noexcept
//      mask representing a mask for a given bit within its word
//
//
//                  TESTS
//
// You can find unit tests for dynamic_bitset in its official repo:
// https://github.com/iboB/itlib/blob/master/test/
//
#pragma once

#include <cstdint>
#include <vector>
#include <cassert>

namespace itlib {

namespace impl {
struct dynamic_bitset_iterator_base
{
    explicit dynamic_bitset_iterator_base(size_t i) noexcept : index(i) {}
    void inc() { ++index; }
    void dec() { --index; }
    bool operator==(const dynamic_bitset_iterator_base& other) const { return index==other.index; }
    bool operator!=(const dynamic_bitset_iterator_base& other) const { return index!=other.index; }
    friend ptrdiff_t operator-(const dynamic_bitset_iterator_base& a, const dynamic_bitset_iterator_base& b) { return ptrdiff_t(a.index) - ptrdiff_t(b.index); }
    size_t index;
};
}

template <typename Buffer = std::vector<uint32_t>>
class dynamic_bitset {
public:
    using buffer_type = Buffer;
    using word_type = typename buffer_type::value_type;
    static_assert(std::is_unsigned<word_type>::value, "word_type must be unsigned");
    using allocator_type = typename buffer_type::allocator_type;
    using size_type = size_t;
    using buffer_size_type = typename buffer_type::size_type;
    static constexpr uint8_t bits_per_word = sizeof(word_type) * 8;

    explicit dynamic_bitset(const allocator_type& alloc = allocator_type()) noexcept : m_buf(alloc), m_size(0) {}
    explicit dynamic_bitset(size_type size, word_type value = 0, const allocator_type& alloc = allocator_type())
        : m_buf(word_size(size), value, alloc), m_size(size)
    {}

    // copy and move
    dynamic_bitset(const dynamic_bitset& x) = default;
    dynamic_bitset& operator=(const dynamic_bitset& x) = default;

    dynamic_bitset(dynamic_bitset&& x) noexcept
        : m_buf(std::move(x.m_buf)), m_size(x.m_size)
    {
        x.m_size = 0;
    }
    dynamic_bitset& operator=(dynamic_bitset&& x) noexcept
    {
        m_buf = std::move(x.m_buf);
        m_size = x.m_size;
        x.m_size = 0;
        return *this;
    }

    // access
    size_type size() const noexcept { return m_size; }
    buffer_size_type word_size() const noexcept { return m_buf.size(); }
    word_type* data() noexcept { return m_buf.data(); }
    const word_type* data() const noexcept { return m_buf.data(); }
    bool empty() const noexcept { return m_buf.empty(); }
    size_t byte_size() const noexcept { return (m_size + 7) / 8; }
    const buffer_type& buffer() const { return m_buf; }

    // test single bit
    bool test(size_type i) const noexcept
    {
        return m_buf[word_index(i)] & word_mask(i);
    }
    // all bits true?
    bool all() const noexcept
    {
        if (empty()) return false;
        auto it = m_buf.begin();
        auto endit = it + m_size / bits_per_word;
        for (; it != endit; ++it)
        {
            if (*it != ~word_type(0)) return false;
        }
        if (endit != m_buf.end())
        {
            word_type mask = word_mask(m_size % bits_per_word) - 1;
            if ((m_buf.back() & mask) != mask) return false;
        }
        return true;
    }
    // any bits true?
    bool any() const noexcept
    {
        if (empty()) return false;
        auto it = m_buf.begin();
        auto endit = it + m_size / bits_per_word;
        for (; it != endit; ++it)
        {
            if (*it != 0) return true;
        }
        if (endit != m_buf.end())
        {
            word_type mask = word_mask(m_size % bits_per_word) - 1;
            if (m_buf.back() & mask ) return true;
        }
        return false;
    }
    // no bits true?
    bool none() const noexcept { return !any(); }

    // modifiers
    // set bit to true
    void set(size_type i)
    {
        m_buf[word_index(i)] |= word_mask(i);
    }
    // set bit to false
    void reset(size_type i)
    {
        m_buf[word_index(i)] &= ~word_mask(i);
    }
    // set bit to value
    void set(size_type i, bool b)
    {
        if (b) set(i);
        else reset(i);
    }
    // flip bit
    void flip(size_type i)
    {
        m_buf[word_index(i)] ^= word_mask(i);
    }

    // size modifiers
    void reserve(size_type size) { m_buf.reserve(word_size(size)); }
    void resize(size_type size)
    {
        auto pad = m_size % bits_per_word;
        if (pad) {
            // clear the tail of the buf with zeroes in case they remain from previous ops
            m_buf.back() &= word_mask(pad) - 1;
        }
        m_size = size;
        m_buf.resize(word_size(size));
    }
    void push_back(bool b)
    {
        resize(m_size + 1);
        set(m_size-1, b);
    }

    void assign(const buffer_type& buf)
    {
        m_buf = buf;
        m_size = m_buf.size() * bits_per_word;
    }

    void append(const dynamic_bitset& other)
    {
        if (other.empty()) return;

        const auto pad = m_size % bits_per_word;
        if (pad == 0)
        {
            // lucky
            m_buf.insert(m_buf.end(), other.m_buf.begin(), other.m_buf.end());
            m_size += other.size();
            return;
        }

        auto back_buf_index = m_buf.size() - 1;
        resize(m_size + other.m_size);
        auto myi = m_buf.begin() + back_buf_index;
        for (size_t i = 0; i<other.m_buf.size(); ++i)
        {
            auto cur = other.m_buf[i];
            *myi |= cur << pad;
            ++myi;
            if (myi == m_buf.end()) break;
            *myi = cur >> (bits_per_word - pad);
        }
    }

    // vector-like
    struct bitref
    {
        bitref(word_type& w, word_type m) : word(w), mask(m) {}
        bool test() const { return word & mask; }
        void set() { word |= mask; }
        void reset() { word &= ~mask; }
        void flip() { word ^= mask; }
        void set(bool b)
        {
            if (b) set();
            else reset();
        }
        operator bool() const { return test(); }
        bitref operator=(bool b) { set(b); return *this; }
    private:
        word_type& word;
        const word_type mask;
    };

    bool operator[](size_type i) const { return test(i); }
    bitref operator[](size_type i) { return bitref(m_buf[word_index(i)], word_mask(i)); }

    struct iterator : public impl::dynamic_bitset_iterator_base
    {
        iterator(dynamic_bitset& b, size_t i) : impl::dynamic_bitset_iterator_base(i), bitset(b) {}
        iterator& operator++() { inc(); return *this; }
        iterator operator++(int) { auto ret = *this; inc(); return ret; }
        iterator& operator--() { dec(); return *this; }
        iterator operator--(int) { auto ret = *this; dec(); return ret; }
        friend iterator operator+(iterator a, ptrdiff_t b) { return iterator(a.bitset, a.index + b); }
        friend iterator operator-(iterator a, ptrdiff_t b) { return iterator(a.bitset, a.index - b); }
        bool operator*() const
        {
            return bitset.test(index);
        }
        bitref operator*()
        {
            return bitset[index];
        }

        dynamic_bitset& bitset;
    };

    struct const_iterator : public impl::dynamic_bitset_iterator_base
    {
        const_iterator(const dynamic_bitset& b, size_t i) : impl::dynamic_bitset_iterator_base(i), bitset(b) {}
        const_iterator(iterator i) : const_iterator(i.bitset, i.index) {}
        const_iterator& operator++() { inc(); return *this; }
        const_iterator operator++(int) { auto ret = *this; inc(); return ret; }
        const_iterator& operator--() { dec(); return *this; }
        const_iterator operator--(int) { auto ret = *this; dec(); return ret; }
        friend const_iterator operator+(const_iterator a, ptrdiff_t b) { return const_iterator(a.bitset, a.index + b); }
        friend const_iterator operator-(const_iterator a, ptrdiff_t b) { return const_iterator(a.bitset, a.index - b); }
        bool operator*() const
        {
            return bitset.test(index);
        }

        const dynamic_bitset& bitset;
    };

    iterator begin() { return iterator(*this, 0); }
    iterator end() { return iterator(*this, m_size); }

    const_iterator begin() const { return const_iterator(*this, 0); }
    const_iterator end() const { return const_iterator(*this, m_size); }

    void copy(iterator from, iterator to, const_iterator source)
    {
        while (from.index % bits_per_word)
        {
            // copy single bits until we reach a word boundary of our bitset
            if (from == to) return;
            *from++ = *source++;
        }

        const auto pad = source.index % bits_per_word;
        auto own_windex = word_index(from.index);
        auto src_windex = word_index(source.index);
        auto& src = source.bitset.m_buf;
        // for each whole word of our bitset
        while (to - from >= bits_per_word)
        {
            m_buf[own_windex] = src[src_windex] >> pad;
            if (pad)
            {
                // there is padding so there must be at least one more element in the source
                m_buf[own_windex] |= src[src_windex + 1] << (bits_per_word - pad);
            }
            ++own_windex;
            ++src_windex;
            from.index += bits_per_word;
        }

        if (from == to) return; // we finished on a word boundary

        // now what remains is some bits from source that need to be splicet onto our last element
        const auto rpad = to.index % bits_per_word;
        m_buf[own_windex] = ((src[src_windex] >> pad) & (word_mask(rpad) - 1)) | (m_buf[own_windex] & ~(word_mask(rpad) - 1));
    }

    // static ops

    // required word size for a given bit size
    static constexpr buffer_size_type word_size(size_type size) noexcept
    {
        return (size + bits_per_word - 1) / bits_per_word;
    }
    // index of a word in buffer for a given bit index
    static constexpr buffer_size_type word_index(size_type i) noexcept
    {
        return i / bits_per_word;
    }
    // index of a bit within its word
    static constexpr uint8_t bit_index(size_type i) noexcept
    {
        return uint8_t(i % bits_per_word);
    }
    // mask representing a mask for a given bit within its word
    static constexpr word_type word_mask(size_type index) noexcept
    {
        return word_type(1) << bit_index(index);
    }

private:
    buffer_type m_buf;
    size_type m_size; // doesn't need to be divisible by bits_per_word
};

}
