// itlib-time_t v1.01
//
// A thin wrapper of std::time_t which provides thread safe std::tm getters and
// type-safe (std::chrono::duration-based) arithmetic
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
//  1.00 (2020-10-36) Initial release
//  1.01 (2021-04-29) Added named ctors: now, from_gmtime, from_localtime
//
//
//                  DOCUMENTATION
//
// Simply include this file wherever you need.
// It defines the class itlib::time_t which is a thin wrapper of std::time_t
//
// It provides multiplatform thread-safe ops to convert to std::tm
// * std::tm itlib::time_t::gmtime() const
// * std::tm itlib::time_t::localtime() const
//
// It also provides type-safe arithmetic
// * itlib::time_t operators +,-,+= and -= with std::chrono_duration
// * std::chrono_duration operator-(itlib::time_t a, itlib::time_t b)
//
// The file also defines the function
// std::string itlib::strftime(const char* fmt, const std::tm& tm);
// It works exactly as std::strftime but returns a std::string with the
// appropriate size
//
//                  TESTS
//
// You can find unit tests for time_t in its official repo:
// https://github.com/iboB/itlib/blob/master/test/
//
#pragma once

#include <ctime>
#include <chrono>
#include <cstdint>
#include <string>

namespace itlib
{

class time_t
{
public:
    using timestamp_type = int64_t;
    using duration_type = std::chrono::duration<timestamp_type>;

    time_t() = default;
    time_t(const time_t&) = default;
    time_t& operator=(const time_t&) = default;

    explicit time_t(const std::time_t& st)
    {
        m_t = static_cast<timestamp_type>(st);
    }
    explicit operator std::time_t() const
    {
        return static_cast<std::time_t>(m_t);
    }

    timestamp_type seconds_since_epoch() const { return m_t; }

    static time_t from_seconds(timestamp_type s) { return time_t(s); }

    static time_t now() { return from_seconds(std::time(nullptr)); }

    // non-const argument - gets normalized internally
    static time_t from_gmtime(std::tm& gmtm)
    {
        auto tt =
#ifdef _WIN32
            _mkgmtime(&gmtm);
#else
            timegm(&gmtm);
#endif
        return from_seconds(tt);
    }

    // non-const argument - gets normalized internally
    static time_t from_localtime(std::tm& localtm)
    {
        return from_seconds(mktime(&localtm));
    }

    // cmp
    friend bool operator==(const time_t& a, const time_t& b) { return a.m_t == b.m_t; }
    friend bool operator!=(const time_t& a, const time_t& b) { return a.m_t != b.m_t; }
    friend bool operator<(const time_t& a, const time_t& b) { return a.m_t < b.m_t; }
    friend bool operator<=(const time_t& a, const time_t& b) { return a.m_t <= b.m_t; }
    friend bool operator>(const time_t& a, const time_t& b) { return a.m_t > b.m_t; }
    friend bool operator>=(const time_t& a, const time_t& b) { return a.m_t >= b.m_t; }

    // arithmetic
    template <typename Rep, typename Period>
    time_t& operator+=(const std::chrono::duration<Rep, Period>& d) { m_t += dc(d); return *this; }
    template <typename Rep, typename Period>
    time_t operator+(const std::chrono::duration<Rep, Period>& d) const { return time_t(m_t + dc(d)); }
    template <typename Rep, typename Period>
    time_t& operator-=(const std::chrono::duration<Rep, Period>& d) { m_t -= dc(d); return *this;  }
    template <typename Rep, typename Period>
    time_t operator-(const std::chrono::duration<Rep, Period>& d) const { return time_t(m_t - dc(d)); }

    friend duration_type operator-(const time_t& a, const time_t& b) { return duration_type(a.m_t - b.m_t); }

    // ops
    std::tm gmtime() const
    {
        std::tm ret = {};
        auto mt = std::time_t(*this);
#if defined(_MSC_VER)
        gmtime_s(&ret, &mt);
#else
        gmtime_r(&mt, &ret);
#endif
        return ret;
    }

    std::tm localtime() const
    {
        std::tm ret = {};
        auto mt = std::time_t(*this);
#if defined(_MSC_VER)
        localtime_s(&ret, &mt);
#else
        localtime_r(&mt, &ret);
#endif
        return ret;
    }

private:
    template <typename Rep, typename Period>
    static timestamp_type dc(const std::chrono::duration<Rep, Period>& d)
    {
        return std::chrono::duration_cast<duration_type>(d).count();
    }

    timestamp_type m_t = 0;
};

inline std::string strftime(const char* format, const std::tm& tm)
{
    std::string ret;
    ret.resize(128);
    size_t len;
    while (!(len = std::strftime(&ret.front(), ret.size(), format, &tm))) ret.resize(2 * ret.size());
    ret.resize(len);
    return ret;
}

}
