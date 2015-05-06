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

#define NOMINMAX
#include "harness_defs.h"
#include "test_concurrent_queue.h"
#include "tbb/concurrent_queue.h"
#include "tbb/concurrent_monitor.cpp"
#include "harness.h"
#include "harness_allocator.h"

#if _MSC_VER==1500 && !__INTEL_COMPILER
    // VS2008/VC9 seems to have an issue; limits pull in math.h
    #pragma warning( push )
    #pragma warning( disable: 4985 )
#endif
#include <limits>
#if _MSC_VER==1500 && !__INTEL_COMPILER
    #pragma warning( pop )
#endif

template <typename Q>
class FloggerBody : NoAssign {
    Q *q;
public:
    FloggerBody(Q *q_) : q(q_) {}
    void operator()(const int threadID) const {
        typedef typename Q::value_type value_type;
        value_type elem = value_type(threadID);
        for (size_t i=0; i<275; ++i) {
            q->push(elem);
            (void) q->try_pop(elem);
        }
    }
};

template <typename HackedQRep, typename Q>
void TestFloggerHelp(HackedQRep* hack_rep, Q* q, size_t items_per_page) {
    size_t nq = HackedQRep::n_queue;
    size_t hack_val = std::numeric_limits<std::size_t>::max() & ~( nq * items_per_page - 1 );
    hack_rep->head_counter = hack_val;
    hack_rep->tail_counter = hack_val;
    size_t k = hack_rep->tail_counter & -(ptrdiff_t)nq;

    for (size_t i=0; i<nq; ++i) {
        hack_rep->array[i].head_counter = k;
        hack_rep->array[i].tail_counter = k;
    }
    NativeParallelFor(MaxThread, FloggerBody<Q>(q));
    ASSERT(q->empty(), "FAILED flogger/empty test.");
    delete q;
}

template <typename T>
void TestFlogger(T /*max*/) {
    {
        tbb::concurrent_queue<T>* q = new tbb::concurrent_queue<T>;
        REMARK("Wraparound on strict_ppl::concurrent_queue...");
        hacked_concurrent_queue_rep* hack_rep = ((hacked_concurrent_queue<T>*)(void*)q)->my_rep;
        TestFloggerHelp(hack_rep, q, hack_rep->items_per_page);
        REMARK(" works.\n");
    }
    {
        tbb::concurrent_bounded_queue<T>* q = new tbb::concurrent_bounded_queue<T>;
        REMARK("Wraparound on tbb::concurrent_bounded_queue...");
        hacked_bounded_concurrent_queue* hack_q = (hacked_bounded_concurrent_queue*)(void*)q;
        hacked_bounded_concurrent_queue_rep* hack_rep = hack_q->my_rep;
        TestFloggerHelp(hack_rep, q, hack_q->items_per_page);
        REMARK(" works.\n");
    }
}

void TestWraparound() {
    REMARK("Testing Wraparound...\n");
    TestFlogger(std::numeric_limits<int>::max());
    TestFlogger(std::numeric_limits<unsigned char>::max());
    REMARK("Done Testing Wraparound.\n");
}

int TestMain () {
    TestWraparound();
    return Harness::Done;
}
