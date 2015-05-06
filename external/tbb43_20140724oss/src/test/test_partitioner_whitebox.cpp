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

#include <algorithm>
#include <typeinfo>

#include "tbb/blocked_range.h"
#include "tbb/tbb_thread.h"
#include "tbb/enumerable_thread_specific.h"

#include "harness.h"
#include "harness_assert.h"
#include "test_partitioner.h"


typedef tbb::enumerable_thread_specific<size_t> ThreadNumsType;
size_t g_threadNumInitialValue = 10;
ThreadNumsType g_threadNums(g_threadNumInitialValue);

// simulate a subset of task.h
namespace tbb {
namespace internal {
typedef unsigned short affinity_id;
}
class fake_task {
public:
    typedef internal::affinity_id affinity_id;
    void set_affinity(affinity_id a) { my_affinity = a; }
    affinity_id affinity() const { return my_affinity; }
    void set_parent(fake_task* p) { my_parent = p; }
    fake_task *parent() const { return my_parent; }
    bool is_stolen_task() const { return false; }
    intptr_t ref_count() const { return 1; }
    bool is_cancelled() const { return false; }
    static void spawn(fake_task &) {} // for legacy in partitioner.h
    virtual fake_task* execute() = 0; // enables dynamic_cast

    fake_task() : my_parent(0), my_affinity(0) {}
    virtual ~fake_task() {}
private:
    fake_task *my_parent;
    affinity_id my_affinity;
};
}

#define __TBB_task_H
#define get_initial_auto_partitioner_divisor my_get_initial_auto_partitioner_divisor
#define affinity_partitioner_base_v3 my_affinity_partitioner_base_v3
#define task fake_task
#define __TBB_STATIC_THRESHOLD 0
#include "tbb/partitioner.h"
#undef __TBB_STATIC_THRESHOLD
#undef task
#undef affinity_partitioner_base_v3
#undef get_initial_auto_partitioner_divisor

#include "string.h"
#include "harness.h"

// replace library functions to simulate concurrency
namespace tbb {
namespace internal {
size_t my_get_initial_auto_partitioner_divisor() {
    const size_t X_FACTOR = 4;
    return X_FACTOR * g_threadNums.local();
}

void* __TBB_EXPORTED_FUNC NFS_Allocate( size_t n_element, size_t element_size, void* hint );
void __TBB_EXPORTED_FUNC NFS_Free( void* );

void my_affinity_partitioner_base_v3::resize( unsigned factor ) {
    // Check factor to avoid asking for number of workers while there might be no arena.
    size_t new_size = factor ? factor * g_threadNums.local() : 0;
    if (new_size != my_size) {
        if (my_array) {
            NFS_Free(my_array);
            // Following two assignments must be done here for sake of exception safety.
            my_array = NULL;
            my_size = 0;
        }
        if (new_size) {
            my_array = static_cast<affinity_id*>(NFS_Allocate(new_size, sizeof(affinity_id), NULL ));
            memset(my_array, 0, sizeof(affinity_id) * new_size);
            my_size = new_size;
        }
    }
}

} //namespace internal
// simulate a subset of parallel_for
namespace interface7 {
namespace internal {

// paralle_for algorithm that executes sequentially
template<typename Range, typename Body, typename Partitioner>
class start_for : public fake_task {
    Range my_range;
    Body my_body;
    typename Partitioner::task_partition_type my_partition;
    size_t m_executedBegin, m_executedEnd;
    bool m_firstTimeRun;
    size_t m_joinedBegin, m_joinedEnd;
    test_partitioner_utils::BinaryTree* m_tree;
public:
    start_for( const Range& range, const Body& body, Partitioner& partitioner,
               test_partitioner_utils::BinaryTree* tree ) :
        my_range(range), my_body(body), my_partition(partitioner),
        m_executedBegin(0), m_executedEnd(0), m_firstTimeRun(true),
        m_joinedBegin(/* grows left */ range.end()), m_joinedEnd(range.end()), m_tree(tree)
    {
        if (m_tree) {
            m_tree->push_node( test_partitioner_utils::make_node(my_range.begin(), my_range.end(), affinity()) );
        }
    }
    //! Splitting constructor used to generate children.
    /** parent_ becomes left child.  Newly constructed object is right child. */
    start_for( start_for& parent_, typename Partitioner::split_type& split_obj) :
        my_range(parent_.my_range, split_obj),
        my_body(parent_.my_body),
        my_partition(parent_.my_partition, split_obj),
        m_executedBegin(0), m_executedEnd(0), m_firstTimeRun(true),
        m_joinedBegin(/* grows left */ my_range.end()), m_joinedEnd(my_range.end()),
        m_tree(parent_.m_tree)
    {
        set_parent(parent_.parent());
        my_partition.set_affinity(*this);

        if (m_tree) {
            // collecting splitting statistics
            m_tree->push_node( test_partitioner_utils::make_node(my_range.begin(), my_range.end(), affinity()) );
            m_tree->push_node( test_partitioner_utils::make_node(parent_.my_range.begin(), parent_.my_range.end(), parent_.affinity()) );
        }
    }
    //! Construct right child from the given range as response to the demand.
    /** parent_ remains left child.  Newly constructed object is right child. */
    start_for( start_for& parent_, const Range& r, depth_t d ) :
        my_range(r),
        my_body(parent_.my_body),
        my_partition(parent_.my_partition, tbb::split()),
        m_executedBegin(0), m_executedEnd(0), m_firstTimeRun(true),
        m_joinedBegin(/* grows left */ r.end()), m_joinedEnd(r.end()),
        m_tree(parent_.m_tree)
    {
        set_parent(parent_.parent());
        my_partition.set_affinity(*this);
        my_partition.align_depth( d );
    }
    fake_task* execute() {
        my_partition.check_being_stolen( *this );
        size_t origBegin = my_range.begin();
        size_t origEnd = my_range.end();

        my_partition.execute(*this, my_range);

        ASSERT(m_executedEnd == m_joinedBegin, "Non-continuous execution");
        m_executedEnd = m_joinedEnd;

        ASSERT(origBegin == m_executedBegin && origEnd == m_executedEnd,
               "Not all iterations were processed");
        return NULL;
    }
    //! Run body for range, serves as callback for partitioner
    void run_body( Range &r ) {
        my_body(r);

        if (m_firstTimeRun) {
            m_firstTimeRun = false;
            m_executedBegin = m_executedEnd = r.begin();
        }

        ASSERT(m_executedBegin <= r.begin() && m_executedEnd < r.end(), "Non-continuous execution");
        m_executedEnd = r.end();
    }
    //! spawn right task, serves as callback for partitioner
    void offer_work(typename Partitioner::split_type& split_obj) {
        start_for sibling(*this, split_obj);
        sibling.execute();
        join(sibling.m_executedBegin, sibling.m_executedEnd);
    }
    //! spawn right task, serves as callback for partitioner
    void offer_work(const Range& r, depth_t d = 0) {
        start_for sibling(*this, r, d);
        sibling.execute();
        join(sibling.m_executedBegin, sibling.m_executedEnd);
    }
    void join(size_t siblingExecutedBegin, size_t siblingExecutedEnd) {
        ASSERT(siblingExecutedEnd == m_joinedBegin, "?");
        m_joinedBegin = siblingExecutedBegin;
    }
};

} //namespace internal
} //namespace interface7
} //namespace tbb

namespace whitebox_simulation {
using namespace tbb::interface7::internal;
template<typename Range, typename Body, typename Partitioner>
void parallel_for( const Range& range, const Body& body, Partitioner& partitioner,
                   test_partitioner_utils::BinaryTree* tree = NULL) {
    if (!range.empty()) {
        flag_task parent;
        start_for<Range, Body, Partitioner> start(range, body, partitioner, tree);
        start.set_parent(&parent);
        start.execute();
    }
}

}//namespace whitebox_simulation

namespace uniform_iterations_distribution {

/*
 * Test checks uniform distribution of range's iterations among all tasks just after
 * work distribution phase has been completed and just before work balancing phase has been started
 */


/*
 * BlockedRange uses tbb::blocked_range formula for proportion calculation
 * InvertedProportionRange inverts proportion suggested by partitioner (e.g. 1:3 --> 3:1)
 * ExactSplitRange uses integer arithmetic for accurate splitting
 */

using namespace test_partitioner_utils;
using tbb::blocked_range;
using tbb::split;
using tbb::proportional_split;

class BlockedRange: public RangeStatisticCollector, public blocked_range<size_t>  {
public:
    BlockedRange(size_t _begin, size_t _end, RangeStatisticData *statData) :
        RangeStatisticCollector(statData),
        blocked_range<size_t>(_begin, _end)
    { }

    BlockedRange(BlockedRange& r, split) :
        RangeStatisticCollector(r, r.size()),
        blocked_range<size_t>(r, split())
    { }

    BlockedRange(BlockedRange& r, proportional_split& p) :
        RangeStatisticCollector(r, p),
        blocked_range<size_t>(r, p)
    { }

    static const bool is_divisible_in_proportion = true;
};

class InvertedProportionRange: public RangeBase<InvertedProportionRange, float> {
public:
    InvertedProportionRange(size_t _begin, size_t _end, RangeStatisticData *statData) :
        RangeBase<InvertedProportionRange, float>(_begin, _end, statData) { }
    InvertedProportionRange(InvertedProportionRange& r, split) :
        RangeBase<InvertedProportionRange, float>(r, split()) { }
    InvertedProportionRange(InvertedProportionRange& r, proportional_split& p) :
        RangeBase<InvertedProportionRange, float>(r, p) { }
    float compute_right_part(RangeBase<InvertedProportionRange, float>& r, tbb::proportional_split& p) {
        return float(r.size() * float(p.left())) / float(p.left() + p.right());
    }
    static const bool is_divisible_in_proportion = true;
};

class ExactSplitRange: public RangeBase<ExactSplitRange, size_t> {
public:
    ExactSplitRange(size_t _begin, size_t _end, RangeStatisticData *statData) :
        RangeBase<ExactSplitRange, size_t>(_begin, _end, statData) { }
    ExactSplitRange(ExactSplitRange& r, split) :
        RangeBase<ExactSplitRange, size_t>(r, split()) { }
    ExactSplitRange(ExactSplitRange& r, proportional_split& p) :
        RangeBase<ExactSplitRange, size_t>(r, p) { }
    size_t compute_right_part(RangeBase<ExactSplitRange, size_t>& r, tbb::proportional_split& p) {
        size_t parts = size_t(p.left() + p.right());
        size_t currSize = r.size();
        size_t int_part = currSize / parts * p.right();
        size_t remainder = currSize % parts * p.right();
        int_part += remainder / parts;
        remainder %= parts;
        size_t right_part = int_part + (remainder > parts/2 ? 1 : 0);
        return right_part;
    }
    static const bool is_divisible_in_proportion = true;
};

template <typename Range, typename Body, typename Partitioner>
void test_case(Range& range, Body& body, Partitioner& partitioner,
               test_partitioner_utils::BinaryTree* tree = NULL) {
    whitebox_simulation::parallel_for(range, body, partitioner, tree);
}

// Functions generate size for range objects used in tests
template <typename T>
size_t default_range_size_generator(T factor, size_t task_num) {
    return size_t(factor * task_num);
}

size_t shifted_left_range_size_generator(size_t factor, size_t task_num) {
    return factor * task_num - 1;
}

size_t shifted_right_range_size_generator(size_t factor, size_t task_num) {
    return factor * task_num + 1;
}

size_t max_range_size_generator(size_t, size_t) {
    return size_t(-1);
}

template <typename RangeType, typename T>
class ParallelTestBody {
public:
    typedef size_t (*RangeSizeGenFunc)(T, size_t);

    ParallelTestBody(T *range_size_factors, unsigned len,
                     size_t tolerance, RangeSizeGenFunc range_size_generator, size_t range_begin) :
        m_factors(range_size_factors), m_len(len), m_tolerance(tolerance),
        m_range_size_generator(range_size_generator), m_range_begin(range_begin)
    { }

    void operator ()(size_t task_num) const {
        task_num++; // since NativeParalleFor starts indexing from zero

        g_threadNums.local() = task_num;

        for (unsigned i = 0; i < m_len; ++i) {
            // initializing
            size_t range_end = m_range_size_generator(m_factors[i], task_num);
            RangeStatisticData stat = {
                0,          // range num
                0,          // minimal size of range
                0,          // maximal size of range
                false       // minimal size of range was not rewritten yet
            };

            RangeType range = RangeType(m_range_begin, range_end, &stat);
            tbb::affinity_partitioner ap;
            SimpleBody body;
            test_case(range, body, ap, NULL);

            // Checking that all threads were given a task
            size_t rangeObjsCreated = stat.m_rangeNum;
            size_t range_size = range_end - m_range_begin;
            if (range_size >= task_num) {
                if (rangeObjsCreated != task_num) {
                    REPORT("During processing of '%s' range of size (%llu) number of range objects created (%lu)"
                           " was not equal to number of tasks (%llu)\n", typeid(range).name(), uint64_t(range_size), uint64_t(rangeObjsCreated), uint64_t(task_num));
                    ASSERT(rangeObjsCreated == task_num, "Incorrect number of range objects was created before work balancing phase started");
                }
            } else if (rangeObjsCreated != range_size && range_size != 0) {
                REPORT("['%s'] number of range objects created (%llu) was not equal to range's size (%llu)"
                       " when number of tasks = %llu\n",
                       typeid(range).name(), uint64_t(rangeObjsCreated), uint64_t(range_size), uint64_t(task_num));
                ASSERT(rangeObjsCreated == range_size, "Incorrect number of range objects was created before work balancing phase started");
            }

            // Checking difference between min and max number of range iterations
            size_t diff = stat.m_maxRangeSize - stat.m_minRangeSize;
            if (diff > m_tolerance) {
                REPORT("Difference (%llu) between maximum and minimum number of splitted '%s' range iterations (%llu), "
                       "while executing test for %llu number of tasks, is greater than allowed tolerance (%llu)\n",
                       uint64_t(diff), typeid(range).name(), uint64_t(range_size), uint64_t(task_num), uint64_t(m_tolerance));
                ASSERT(diff <= m_tolerance, "Uniform iteration distribution error");
            }
        }
    }
private:
    T *m_factors;                              // array of multipliers for parallel_for's number of iterations
    unsigned m_len;                            // array size
    size_t m_tolerance;                        // tolerance for min and max size of ranges
    RangeSizeGenFunc m_range_size_generator;   // function generates range sizes
    size_t m_range_begin;                      // beginning of range iterations
};

template <typename RangeType, typename T>
void test(T range_size_factors[], unsigned len, size_t tolerance = 2,
          size_t (*rsgFunc)(T, size_t) = default_range_size_generator<T>, size_t range_begin = 0) {
    // Ideally difference should be equal to zero in case if iterations can be divided
    // without a remainder and equals to one otherwise but due to floating point rounding
    // on some systems difference reaches default value of 'tolerance' arg

    // some reasonable value for possible number of threads
    size_t max_simulated_threads = 1024;
    size_t hw_threads_num = tbb::tbb_thread::hardware_concurrency();
    size_t threadsToRunOn = std::min<size_t>(max_simulated_threads, hw_threads_num);

    for (size_t task_num = 1; task_num <= max_simulated_threads; task_num += threadsToRunOn) {
        ParallelTestBody<RangeType, T> b = ParallelTestBody<RangeType, T>
            (range_size_factors, len, tolerance, rsgFunc, range_begin);
        NativeParallelFor(threadsToRunOn, b);
    }
}

void test() {
    using namespace test_partitioner_utils::TestRanges;

    {
        // Multipliers for number of tasks (iterations should be distributed uniformly)
        size_t range_size_factors[] = { 1, 2, 3, 4, 5, 7, 9, 13, 27, 29, 30, 31, 32 };
        unsigned len = sizeof(range_size_factors) / sizeof(range_size_factors[0]);

        // tolerance for blocked_range is equal to zero
        test<BlockedRange>(range_size_factors, len, 0);
        test<InvertedProportionRange>(range_size_factors, len);
        test<RoundedDownRange>(range_size_factors, len);
        test<RoundedUpRange>(range_size_factors, len);

        // check only quantity of range objects
        test<Range1_2>(range_size_factors, len, size_t(-1));
        test<Range1_999>(range_size_factors, len, size_t(-1));
        test<Range999_1>(range_size_factors, len, size_t(-1));
    }

    {
        // Multipliers for number of tasks (iterations might not be distributed uniformly)
        float range_size_factors[] = { 1.2f, 2.5f, 3.7f, 4.2f, 5.1f, 8.9f, 27.8f };
        unsigned len = sizeof(range_size_factors) / sizeof(range_size_factors[0]);

        // tolerance for blocked_range is equal to one
        test<BlockedRange>(range_size_factors, len, 1);
        test<InvertedProportionRange>(range_size_factors, len);
        test<RoundedDownRange>(range_size_factors, len);
        test<RoundedUpRange>(range_size_factors, len);
    }


    {
        // Multipliers for number of tasks (iterations might not be distributed uniformly)
        size_t range_size_factors[] = { 1, 2, 3, 4, 5, 7, 9, 11, 13, 27, 29, 30, 31, 32 };
        unsigned len = sizeof(range_size_factors) / sizeof(range_size_factors[0]);

        // tolerance for blocked_range is equal to one
        test<BlockedRange>(range_size_factors, len, 1, &shifted_left_range_size_generator);
        test<BlockedRange>(range_size_factors, len, 1, &shifted_right_range_size_generator);

        test<InvertedProportionRange>(range_size_factors, len, 2, &shifted_left_range_size_generator);
        test<InvertedProportionRange>(range_size_factors, len, 2, &shifted_right_range_size_generator);

        test<RoundedDownRange>(range_size_factors, len, 2, &shifted_left_range_size_generator);
        test<RoundedDownRange>(range_size_factors, len, 2, &shifted_right_range_size_generator);

        test<RoundedUpRange>(range_size_factors, len, 2, &shifted_left_range_size_generator);
        test<RoundedUpRange>(range_size_factors, len, 2, &shifted_right_range_size_generator);
    }


    {
        size_t range_size_factors[] = { 0 };
        unsigned len = sizeof(range_size_factors) / sizeof(range_size_factors[0]);

        // tolerance is 1 since range iterations number is not divided without a remainder
        test<ExactSplitRange>(range_size_factors, len, 1, &max_range_size_generator);
        test<ExactSplitRange>(range_size_factors, len, 1, &max_range_size_generator, size_t(-1) - 10000);
    }
}

} // namespace uniform_iterations_distribution

namespace overflow_during_split {

using tbb::blocked_range;
using tbb::proportional_split;
using tbb::interface7::internal::adaptive_partition_type_base;
using tbb::internal::uint64_t;

class partitioner: public adaptive_partition_type_base<partitioner> {
    size_t m_right_part;
public:
    partitioner(size_t rp) : adaptive_partition_type_base<partitioner>(), m_right_part(rp) { }
    partitioner(partitioner &src, const proportional_split &p) :
        adaptive_partition_type_base<partitioner>(src, p) {
        // computation without potential overflow
        size_t parts = p.left() + p.right();
        size_t int_part = g_threadNums.local() / parts;
        size_t fraction = g_threadNums.local() - int_part * parts; // fraction < parts
        size_t right_divisor = int_part * src.m_right_part + fraction * src.m_right_part / parts;

        // Division in 'right_divisor' very likely is inexact also.
        size_t tolerance = 1;
        size_t diff = (right_divisor < my_divisor) ? (my_divisor - right_divisor) : (right_divisor - my_divisor);
        if (diff > tolerance) {
            REPORT("Difference between %llu and %llu is >= 2, but should be <\n", uint64_t(my_divisor),
                   uint64_t(right_divisor));
            ASSERT(diff <= tolerance, "Overflow occurred in 'adaptive_partition_type_base'");
        }
    }
};

void test() {
    g_threadNums.local() = size_t(-1) / 4;
    size_t right_part = 6;
    partitioner fp(right_part);
    // trying to generate overflow in adaptive_partition_type_base
    partitioner(fp, proportional_split(2, right_part));
}

} // namespace overflow_during_split

// In order to see splitting in action enable 'VISUALIZE_SPLITTING'
// macro alogn with the following options:
// - PARTITIONER='partitioner-type' - partitioner to be used for splitting
// - RANGE_SIZE='<number>' - number of initial iterations in tbb::blocked_range
// - THREADS_NUM='<number>' - simulated number of threads to be used
// If these options are not defined the default value will be used:
// PARTITIONER = tbb::affinity_partitioner, RANGE_SIZE = 123, THREADS_NUM = 10
void visualize_splitting_topology() {
    using test_partitioner_utils::BinaryTree;
    using test_partitioner_utils::SimpleBody;
    using uniform_iterations_distribution::test_case;
    using tbb::blocked_range;

#ifndef PARTITIONER
#define PARTITIONER tbb::affinity_partitioner
#endif
#define PARTITIONER_NAME2(p) #p
#define PARTITIONER_NAME(p) PARTITIONER_NAME2(p)
#ifndef RANGE_SIZE
#define RANGE_SIZE 123
#endif
#ifndef THREADS_NUM
#define THREADS_NUM 10
#endif

    g_threadNums.local() = THREADS_NUM;
    PARTITIONER partitioner;
    blocked_range<size_t> range(0, RANGE_SIZE);
    BinaryTree t;
    SimpleBody body;
    test_case(range, body, partitioner, &t);
    REPORT("partitioner: %s\nrange_size: %u\nthreads num: %u\n", PARTITIONER_NAME(PARTITIONER), unsigned(RANGE_SIZE), unsigned(THREADS_NUM));
    t.visualize();
}

int TestMain () {
#ifdef VISUALIZE_SPLITTING
    visualize_splitting_topology();
#else
    uniform_iterations_distribution::test();
    overflow_during_split::test();
#endif

    return Harness::Done;
}
