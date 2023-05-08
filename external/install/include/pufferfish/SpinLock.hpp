//
// RapMap - Rapid and accurate mapping of short reads to transcriptomes using
// quasi-mapping.
// Copyright (C) 2015, 2016 Rob Patro, Avi Srivastava, Hirak Sarkar
//
// This file is part of RapMap.
//
// RapMap is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// RapMap is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RapMap.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef __SPIN_LOCK_HPP__
#define __SPIN_LOCK_HPP__

#include <atomic>

// Taken from http://stackoverflow.com/questions/26583433/c11-implementation-of-spinlock-using-atomic
class SpinLock {
    std::atomic_flag locked = ATOMIC_FLAG_INIT ;
public:
    void lock() {
        while (locked.test_and_set(std::memory_order_acquire)) { ; }
    }

    // from http://stackoverflow.com/questions/19742993/implementing-a-spinlock-in-boost-example-neededhttp://stackoverflow.com/questions/19742993/implementing-a-spinlock-in-boost-example-needed
    // is this legit?
    bool try_lock() {
        return !locked.test_and_set(std::memory_order_acquire);
    }

    void unlock() {
        locked.clear(std::memory_order_release);
    }
};

#endif //__SPIN_LOCK_HPP__
