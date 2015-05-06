/*
    Copyright 2005-2014 Intel Corporation.  All Rights Reserved.

    This file is part of Threading Building Blocks. Threading Building Blocks is free software;
    you can redistribute it and/or modify it under the terms of the GNU General Public License
    version 2  as  published  by  the  Free Software Foundation.  Threading Building Blocks is
    distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
    implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    See  the GNU General Public License for more details.   You should have received a copy of
    the  GNU General Public License along with Threading Building Blocks; if not, write to the
    Free Software Foundation, Inc.,  51 Franklin St,  Fifth Floor,  Boston,  MA 02110-1301 USA

    As a special exception,  you may use this file  as part of a free software library without
    restriction.  Specifically,  if other files instantiate templates  or use macros or inline
    functions from this file, or you compile this file and link it with other files to produce
    an executable,  this file does not by itself cause the resulting executable to be covered
    by the GNU General Public License. This exception does not however invalidate any other
    reasons why the executable file might be covered by the GNU General Public License.
*/

#include "tbb/atomic.h"
#include "tbb/cache_aligned_allocator.h"
#include "tbb/spin_mutex.h"
#include "tbb/tbb_machine.h"

#include "tbb/concurrent_monitor.h"

struct hacked_micro_queue {
    tbb::atomic<uintptr_t> head_page;
    tbb::atomic<size_t> head_counter;

    tbb::atomic<uintptr_t> tail_page;
    tbb::atomic<size_t> tail_counter;

    tbb::spin_mutex page_mutex;
 };

// hacks for strict_ppl::concurrent_queue
struct hacked_concurrent_queue_rep {
    static const size_t phi = 3;
    static const size_t n_queue = 8;

    tbb::atomic<size_t> head_counter;
    char pad1[tbb::internal::NFS_MaxLineSize-sizeof(tbb::atomic<size_t>)];
    tbb::atomic<size_t> tail_counter;
    char pad2[tbb::internal::NFS_MaxLineSize-sizeof(tbb::atomic<size_t>)];

    size_t items_per_page;
    size_t item_size;
    tbb::atomic<size_t> n_invalid_entries;
    char pad3[tbb::internal::NFS_MaxLineSize-sizeof(size_t)-sizeof(size_t)-sizeof(tbb::atomic<size_t>)];

    hacked_micro_queue array[n_queue];
};

struct hacked_concurrent_queue_page_allocator {
    size_t foo;
};

template <typename T>
struct hacked_concurrent_queue : public hacked_concurrent_queue_page_allocator {
    hacked_concurrent_queue_rep* my_rep;
    typename tbb::cache_aligned_allocator<T>::template rebind<char>::other my_allocator;
};

// hacks for concurrent_bounded_queue and deprecated::concurrent_queue
struct hacked_bounded_concurrent_queue_rep {
    static const size_t phi = 3;
    static const size_t n_queue = 8;

    tbb::atomic<size_t> head_counter;
    char cmon_items_avail[ sizeof(tbb::internal::concurrent_monitor) ];
    tbb::atomic<size_t> n_invalid_entries;
    char pad1[tbb::internal::NFS_MaxLineSize-((sizeof(tbb::atomic<size_t>)+sizeof(tbb::internal::concurrent_monitor)+sizeof(tbb::atomic<size_t>))&(tbb::internal::NFS_MaxLineSize-1))];

    tbb::atomic<size_t> tail_counter;
    char cmon_slots_avail[ sizeof(tbb::internal::concurrent_monitor) ];
    char pad2[tbb::internal::NFS_MaxLineSize-((sizeof(tbb::atomic<size_t>)+sizeof(tbb::internal::concurrent_monitor))&(tbb::internal::NFS_MaxLineSize-1))];
    hacked_micro_queue array[n_queue];

    static const ptrdiff_t infinite_capacity = ptrdiff_t(~size_t(0)/2);
};

struct hacked_bounded_concurrent_queue  {
    size_t foo;
    hacked_bounded_concurrent_queue_rep* my_rep;
    ptrdiff_t my_capacity;
    size_t items_per_page;
    size_t item_size;
};
