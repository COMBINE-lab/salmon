// itlib-flat-set v1.02
//
// std::set-like class with an underlying vector
//
// SPDX-License-Identifier: MIT
// MIT License:
// Copyright(c) 2021 Borislav Stanimirov
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
//  1.02 (2021-09-28) Fixed construction from std::initializer_list which
//                    allowed duplicate elements to find their wey in the set
//  1.01 (2021-09-15) Constructors from std::initializer_list
//  1.00 (2021-08-10) Initial-release
//
//
//                  DOCUMENTATION
//
// Simply include this file wherever you need.
// It defines the class itlib::flat_set, which is an almsot drop-in replacement
// of std::set. Flat set has an optional underlying container which by default
// is std::vector. Thus the items in the set are in a continuous block of
// memory. Thus iterating over the set is cache friendly, at the cost of
// O(n) for insert and erase.
//
// The elements inside (like in std::set) are kept in an order sorted by key.
// Getting a value by key is O(log2 n)
//
// It generally performs much faster than std::set for smaller sets of elements
//
// The difference with std::set, which makes flat_set an not-exactly-drop-in
// replacement is the last template argument:
// * std::set has <key, compare, allocator>
// * itlib::flat_set has <key, compare, container>
// The container must be an std::vector compatible type (itlib::static_vector
// is, for example, viable). The container value type must be `key`
//
//                  Changing the allocator.
//
// If you want to change the allocator of flat set, you'll have to provide a
// container with the appropriate one. Example:
//
// itlib::flat_set<
//      string,
//      less<string>,
//      std::vector<string>, MyAllocator<string>>
//  > myset
//
//
//                  Configuration
//
// itlib::flat_set has a single configurable setting:
//
// const char* overloads
// By default itlib::flat_set provides overloads for the access methods
// (find, lower_bound, count) for const char* for cases when
// std::string is the key, so that no allocations happen when accessing with
// a C-string of a string literal.
// However if const char* or any other class with implicit conversion from
// const char* is the key, they won't compile.
// If you plan on using flat_set with such keys, you'll need to define
// ITLIB_FLAT_SET_NO_CONST_CHAR_OVERLOADS before including the header
//
//
//                  TESTS
//
// You can find unit tests for static_vector in its official repo:
// https://github.com/iboB/itlib/blob/master/test/
//
#pragma once

#include <vector>
#include <algorithm>
#include <type_traits>

#if !defined(ITLIB_FLAT_SET_NO_CONST_CHAR_OVERLOADS)
#include <cstring>
#endif

namespace itlib
{

template <typename Key, typename Compare = std::less<Key>, typename Container = std::vector<Key>>
class flat_set
{
public:
    using key_type = Key;
    using value_type = Key;
    using container_type = Container;
    using key_compare = Compare;
    using reference = value_type&;
    using const_reference = const value_type& ;
    using allocator_type = typename container_type::allocator_type;
    using pointer = typename std::allocator_traits<allocator_type>::pointer;
    using const_pointer = typename std::allocator_traits<allocator_type>::pointer;
    using iterator = typename container_type::iterator;
    using const_iterator = typename container_type::const_iterator;
    using reverse_iterator = typename container_type::reverse_iterator;
    using const_reverse_iterator = typename container_type::const_reverse_iterator;
    using difference_type = typename container_type::difference_type;
    using size_type = typename container_type::size_type;

    flat_set()
    {}

    explicit flat_set(const key_compare& comp, const allocator_type& alloc = allocator_type())
        : m_cmp(comp)
        , m_container(alloc)
    {}

    explicit flat_set(container_type container, const key_compare& comp = key_compare())
        : m_cmp(comp)
        , m_container(std::move(container))
    {
        std::sort(m_container.begin(), m_container.end(), m_cmp);
        auto new_end = std::unique(m_container.begin(), m_container.end());
        m_container.erase(new_end, m_container.end());
    }

    flat_set(std::initializer_list<value_type> init, const key_compare& comp = key_compare(), const allocator_type& alloc = allocator_type())
        : flat_set(container_type(std::move(init), alloc), comp)
    {}

    flat_set(std::initializer_list<value_type> init, const allocator_type& alloc)
        : flat_set(std::move(init), key_compare(), alloc)
    {}

    flat_set(const flat_set& x) = default;
    flat_set(flat_set&& x) = default;

    flat_set& operator=(const flat_set& x) = default;
    flat_set& operator=(flat_set&& x) = default;

    iterator begin() noexcept { return m_container.begin(); }
    const_iterator begin() const noexcept { return m_container.begin(); }
    iterator end() noexcept { return m_container.end(); }
    const_iterator end() const noexcept { return m_container.end(); }
    reverse_iterator rbegin() noexcept { return m_container.rbegin(); }
    const_reverse_iterator rbegin() const noexcept { return m_container.rbegin(); }
    reverse_iterator rend() noexcept { return m_container.rend(); }
    const_reverse_iterator rend() const noexcept { return m_container.rend(); }
    const_iterator cbegin() const noexcept { return m_container.cbegin(); }
    const_iterator cend() const noexcept { return m_container.cend(); }

    bool empty() const noexcept { return m_container.empty(); }
    size_type size() const noexcept { return m_container.size(); }
    size_type max_size() const noexcept { return m_container.max_size(); }

    void reserve(size_type count) { return m_container.reserve(count); }
    size_type capacity() const noexcept { return m_container.capacity(); }

    void clear() noexcept { m_container.clear(); }

    iterator lower_bound(const key_type& k)
    {
        return std::lower_bound(m_container.begin(), m_container.end(), k, m_cmp);
    }

    const_iterator lower_bound(const key_type& k) const
    {
        return std::lower_bound(m_container.begin(), m_container.end(), k, m_cmp);
    }

    iterator find(const key_type& k)
    {
        auto i = lower_bound(k);
        if (i != end() && !m_cmp(k, *i))
            return i;

        return end();
    }

    const_iterator find(const key_type& k) const
    {
        auto i = lower_bound(k);
        if (i != end() && !m_cmp(k, *i))
            return i;

        return end();
    }

    size_t count(const key_type& k) const
    {
        return find(k) == end() ? 0 : 1;
    }

    template <typename P>
    std::pair<iterator, bool> insert(P&& val)
    {
        auto i = lower_bound(val);
        if (i != end() && !m_cmp(val, *i))
        {
            return { i, false };
        }

        return{ m_container.emplace(i, std::forward<P>(val)), true };
    }

    std::pair<iterator, bool> insert(const value_type& val)
    {
        auto i = lower_bound(val);
        if (i != end() && !m_cmp(val, *i))
        {
            return { i, false };
        }

        return{ m_container.emplace(i, val), true };
    }

    template <typename... Args>
    std::pair<iterator, bool> emplace(Args&&... args)
    {
        value_type val(std::forward<Args>(args)...);
        return insert(std::move(val));
    }

    iterator erase(const_iterator pos)
    {
        return m_container.erase(pos);
    }

    size_type erase(const key_type& k)
    {
        auto i = find(k);
        if (i == end())
        {
            return 0;
        }

        erase(i);
        return 1;
    }

    void swap(flat_set& x)
    {
        std::swap(m_cmp, x.m_cmp);
        m_container.swap(x.m_container);
    }

    const container_type& container() const noexcept
    {
        return m_container;
    }

    // DANGER! If you're not careful with this function, you may irreversably break the set
    container_type& modify_container() noexcept
    {
        return m_container;
    }

#if !defined(ITLIB_FLAT_SET_NO_CONST_CHAR_OVERLOADS)
    ///////////////////////////////////////////////////////////////////////////////////
    // const char* overloads for sets with an std::string key to avoid allocs
    iterator lower_bound(const char* k)
    {
        static_assert(std::is_same<std::string, key_type>::value, "flat_set::lower_bound(const char*) works only for std::strings");
        static_assert(std::is_same<std::less<std::string>, key_compare>::value, "flat_set::lower_bound(const char*) works only for std::string-s, compared with std::less<std::string>");
        return std::lower_bound(m_container.begin(), m_container.end(), k, [](const value_type& a, const char* b) -> bool
        {
            return strcmp(a.c_str(), b) < 0;
        });
    }

    const_iterator lower_bound(const char* k) const
    {
        static_assert(std::is_same<std::string, key_type>::value, "flat_set::lower_bound(const char*) works only for std::strings");
        static_assert(std::is_same<std::less<std::string>, key_compare>::value, "flat_set::lower_bound(const char*) works only for std::string-s, compared with std::less<std::string>");
        return std::lower_bound(m_container.begin(), m_container.end(), k, [](const value_type& a, const char* b) -> bool
        {
            return strcmp(a.c_str(), b) < 0;
        });
    }

    iterator find(const char* k)
    {
        auto i = lower_bound(k);
        if (i != end() && *i == k)
            return i;

        return end();
    }

    const_iterator find(const char* k) const
    {
        auto i = lower_bound(k);
        if (i != end() && *i == k)
            return i;

        return end();
    }

    size_t count(const char* k) const
    {
        return find(k) == end() ? 0 : 1;
    }

#endif // !defined(ITLIB_FLAT_SET_NO_CONST_CHAR_OVERLOADS)

private:
    key_compare m_cmp;
    container_type m_container;
};

template <typename Key, typename Compare, typename Container>
bool operator==(const flat_set<Key, Compare, Container>& a, const flat_set<Key, Compare, Container>& b)
{
    return a.container() == b.container();
}

template <typename Key, typename Compare, typename Container>
bool operator!=(const flat_set<Key, Compare, Container>& a, const flat_set<Key, Compare, Container>& b)
{
    return a.container() != b.container();
}

}
