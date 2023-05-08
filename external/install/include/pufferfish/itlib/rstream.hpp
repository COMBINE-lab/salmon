// itlib-rstream v1.00
//
// std::stream-like classes which impose more restrictions on reading
// thus allowing somewhat optimal stream redirection
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
//  1.00 (2020-10-21) Initial release
//
//
//                  DOCUMENTATION
//
// Simply include this file wherever you need.
// It defines two classes
// itlib::rstream (read stream) which wraps a std::istream and allows reads,
// but does not allow seeks. Thus can be sure that all ops in the stream are
// sequential and don't have to worry that seeks might mess-up your reades.
//
// This allows us to define the second class which is called redirect_rstream
// it allows you to to have a rstream from a specific pos within an existing
// std::stream. redirect_rstream seeks to the position when created and seeks
// back to the original position within the wrapped std::stream when destroyed
//
// Configuration:
// You can optionally define ITLIB_RSTREAM_OVERLOAD_RSHIFT to make rstream
// provide a `>>` operator like the one in std::istream.
//
// It is not enabled by default, because it is DANGEROUS when used with
// redirect_stream
//
// The point of having a redirect_stream is to potentially wrap several objects
// which you can read from an rstream into a single block (file, memory...)
// Then, if the packing is tight, reading the last element with >> from an
// object can overflow into the next. Only allow `operator >>` if you know what
// you're doing.
//
//
//                  TESTS
//
// You can find unit tests for rstream in its official repo:
// https://github.com/iboB/itlib/blob/master/test/
//
#pragma once

#include <istream>

namespace itlib
{

template <typename Stream>
class basic_rstream
{
public:
    using stream_type = Stream;
    using char_type = typename Stream::char_type;

    basic_rstream(Stream& in)
        : m_in(in)
    {}

    virtual ~basic_rstream() = default;

    basic_rstream& read(char_type* buf, std::streamsize count)
    {
        m_in.read(buf, count);
        return *this;
    }

    bool fail() const { return m_in.fail(); }
    explicit operator bool() const { return !fail(); }
    bool good() const { return m_in.good(); }
    bool eof() const { return m_in.eof(); }

#if ITLIB_RSTREAM_OVERLOAD_RSHIFT
    template <typename T>
    friend basic_rstream& operator>>(basic_rstream& in, T& t)
    {
        in.m_in >> t;
        return in;
    }
#endif

protected:
    Stream& m_in;
};

using rstream = basic_rstream<std::istream>;

template <typename Stream>
class basic_redirect_rstream : public basic_rstream<Stream>
{
public:
    basic_redirect_rstream(Stream& in, std::streampos pos)
        : basic_rstream<Stream>(in)
        , m_restore(pos)
    {
        m_restore = in.tellg();
        in.seekg(pos);
    }

    ~basic_redirect_rstream()
    {
        this->m_in.seekg(m_restore);
    }
private:
    std::streampos m_restore;
};

using redirect_rstream = basic_redirect_rstream<std::istream>;

}
