// itlib-mem-streambuf v1.01
//
// std::streambuf implementations for working with contiguous memory
//
// SPDX-License-Identifier: MIT
// MIT License:
// Copyright(c) 2020 Borislav Stanimirov
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
//  1.01 (2021-11-18) Fixed mem_ostreambuf bug when used with containers whose
//                    data() returns non-null when empty
//  1.00 (2020-10-16) Initial release
//
//
//                  DOCUMENTATION
//
// Simply include this file wherever you need.
// It defines two classes which help working with std::stream-s with contguous
// memory.
//
// *** itlib::mem_ostreambuf ***
//
// Works with std::ostream. Has a template container which hold contiguous
// memory. The container must provide methods reserve, resize and data.
// Suitable containers are std::vector, std::string and itlib's small_vector
// or static_vector. It will grow the container appropriately so it holds the
// data which was written to the stream.
//
// mem_ostreambuf::peek_container
// Allows you to peek at the current data inside of the container
//
// mem_ostreambuf::get_container
// Transfer ownership of the container to the caller and clears the buffer
//
// Example:
//
// itlib::mem_ostreambuf<std::string> buf;
// std::ostream out(&buf);
// out << "Hello world!";
// auto str = buf.get_container();
// assert(str == "Hello world!");
//
// *** itlib::mem_istreambuf ***
//
// Works with std::istream. Works with a buffer of a given size provided by
// pointer and size. It makes no allocations whatsoever and doesn't touch
// the provided buffer.
//
// Example:
//
// std::string input = "1 2 3";
// itlib::mem_istreambuf<char> buf(input.data(), input.length());
// std::istream in(&buf);
// int a, b, c;
// in >> a >> b >> c;
// assert(a == 1);
// assert(b == 2);
// assert(c == 3);
//
//
//                  TESTS
//
// You can find unit tests for mem-streambuf in its official repo:
// https://github.com/iboB/itlib/blob/master/test/
//
#pragma once

#include <streambuf>
#include <cstring>
#include <cassert>

namespace itlib
{

template <typename Container>
class mem_ostreambuf final : public std::basic_streambuf<typename Container::value_type>
{
private:
    Container m_data;
    static_assert(std::is_pod<typename Container::value_type>::value, "mem ostream must be of pod type");
    using super = std::basic_streambuf<typename Container::value_type>;
public:
    using int_type = typename super::int_type;
    using char_type = typename super::char_type;
    using pos_type = typename super::pos_type;
    using off_type = typename super::off_type;

    mem_ostreambuf(size_t reserve = 0)
    {
        this->setp(m_data.data(), m_data.data());
        cap_resize_by(reserve);
    }

    // put offset
    size_t poff() const
    {
        return size_t(this->pptr() - m_data.data());
    }

    const Container& peek_container() const
    {
        return m_data;
    }

    Container get_container()
    {
        Container ret;
        auto size = this->poff();
        ret.swap(m_data);
        this->setp(nullptr, nullptr);
        ret.resize(size);
        return ret;
    }

    void clear()
    {
        m_data.clear();
        this->setp(m_data.data(), m_data.data());
    }

private:

    int_type overflow(int_type ch) override
    {
        cap_resize_by(1);

        *this->pptr() = char_type(ch);
        this->pbump(1);

        return ch;
    }

    std::streamsize xsputn(const char_type* s, std::streamsize num) override
    {
        // hacky check of size
        if (this->pptr() + num > this->epptr())
        {
            cap_resize_by(this->pptr() + num - this->epptr());
        }

        memcpy(this->pptr(), s, num * sizeof(char_type));
        this->pbump(int_type(num));

        return num;
    }

    pos_type seekpos(pos_type sp, std::ios_base::openmode) override
    {
        this->setp(m_data.data() + int(sp), m_data.data() + m_data.size());
        return sp;
    }

    pos_type seekoff(off_type off, std::ios_base::seekdir way, std::ios_base::openmode) override
    {
        if (way == std::ios_base::cur) return seekpos(poff() + off, std::ios_base::out);
        if (way == std::ios_base::beg) return seekpos(off, std::ios_base::out);
        if (way == std::ios_base::end) return seekpos(m_data.size() + off, std::ios_base::out);
        return super::traits_type::eof();
    }

private:

    void cap_resize_by(size_t by)
    {
        auto off = poff();

        // we need two resize calls to adopt the container growth factor
        m_data.resize(m_data.size() + by);
        auto new_size = m_data.capacity();
        m_data.resize(new_size);

        this->setp(m_data.data() + off, m_data.data() + new_size);
    }
};

template <typename CharT>
class mem_istreambuf final : public std::basic_streambuf<CharT>
{
private:
    static_assert(std::is_pod<CharT>::value, "mem ostream must be of pod type");
    using super = std::basic_streambuf<CharT>;
public:
    using int_type = typename super::int_type;
    using char_type = typename super::char_type;
    using pos_type = typename super::pos_type;
    using off_type = typename super::off_type;

    mem_istreambuf() = default;

    mem_istreambuf(const CharT* beg, size_t size)
    {
        reset(beg, size);
    }

    void reset(const CharT* beg, size_t size)
    {
        auto ucbeg = const_cast<CharT*>(beg);
        this->setg(ucbeg, ucbeg, ucbeg + size);
    }

    // get offset
    size_t goff() const
    {
        return this->gptr() - this->eback();
    }

    size_t size() const
    {
        return this->egptr() - this->eback();
    }

private:
    pos_type seekpos(pos_type sp, std::ios_base::openmode) override
    {
        this->setg(this->eback(), this->eback() + sp, this->egptr());
        return sp;
    }

    pos_type seekoff(off_type off, std::ios_base::seekdir way, std::ios_base::openmode) override
    {
        if (way == std::ios_base::cur) return seekpos(goff() + off, std::ios_base::in);
        if (way == std::ios_base::beg) return seekpos(off, std::ios_base::in);
        if (way == std::ios_base::end) return seekpos(size() + off, std::ios_base::in);
        return super::traits_type::eof();
    }

    std::streamsize xsgetn(char_type* s, std::streamsize count) override
    {
        if (this->gptr() + count > this->egptr())
        {
            count = this->egptr() - this->gptr();
        }

        memcpy(s, this->gptr(), count * sizeof(char_type));
        this->setg(this->eback(), this->gptr() + count, this->egptr());

        return count;
    }
};

}
