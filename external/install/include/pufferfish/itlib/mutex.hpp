// itlib-mutex v0.01 alpha
//
// Mutex types to extend the existing standard mutexes
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
//  0.01 (2021-10-24) Initial release
//
//
//                  DOCUMENTATION
//
// Simply include this file wherever you need.
// It defines the class several mutex types to extend the existing standard
// mutexes with additional features
//
// Nomenclature:
//
// The types are named by concatenating features with underscore with "mutex"
// at the end. The features are:
// * try: blocking is not supported. Only try_lock
// * strong: no spurious behavior. try_lock guarantees a lock exists
// * rec: recursive. Recursive locks from the same thread are supported
//
// Defined types:
// * strong_try_rec_mutex
//   Basically the only difference between this and std::recursive_mutex is the
//   guarantee that try_lock has no spurious false returns
//
//                  TESTS
//
// You can find unit tests for ufunction in its official repo:
// https://github.com/iboB/itlib/blob/master/test/
//
#pragma once

#include <thread>
#include <mutex>
#include <cassert>

namespace itlib
{

class strong_try_rec_mutex
{
public:
    bool try_lock()
    {
        std::lock_guard<std::mutex> l(m_mutex);
        auto tid = std::this_thread::get_id();
        if (m_owner == std::thread::id())
        {
            assert(m_depth == 0);
            m_owner = tid;
        }
        else if (tid != m_owner)
        {
            return false;
        }
        ++m_depth;
        return true;
    }

    void unlock()
    {
        std::lock_guard<std::mutex> l(m_mutex);
        assert(m_owner == std::this_thread::get_id());
        assert(m_depth > 0);
        --m_depth;
        if (m_depth == 0)
        {
            m_owner = std::thread::id();
        }
    }


private:
    std::mutex m_mutex;
    std::thread::id m_owner = std::thread::id(); // current owner
    std::int_fast32_t m_depth = 0; // recursion depth
};

}
