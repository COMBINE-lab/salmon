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

#include "tbb/enumerable_thread_specific.h"
#include "tbb/task_scheduler_init.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_reduce.h"
#include "tbb/blocked_range.h"
#include "tbb/tick_count.h"
#include "tbb/tbb_allocator.h"
#include "tbb/tbb_thread.h"
#include "tbb/atomic.h"

#if !TBB_USE_EXCEPTIONS && _MSC_VER
    // Suppress "C++ exception handler used, but unwind semantics are not enabled" warning in STL headers
    #pragma warning (push)
    #pragma warning (disable: 4530)
#endif

#include <cstring>
#include <vector>
#include <deque>
#include <list>
#include <map>
#include <utility>

#if !TBB_USE_EXCEPTIONS && _MSC_VER
    #pragma warning (pop)
#endif

#include "harness_assert.h"
#include "harness.h"

#if __TBB_GCC_WARNING_SUPPRESSION_PRESENT
#pragma GCC diagnostic ignored "-Wuninitialized"
#endif

static tbb::atomic<int> construction_counter;
static tbb::atomic<int> destruction_counter;

#if TBB_USE_DEBUG
const int REPETITIONS = 4;
const int N = 10000;
const int RANGE_MIN=1000;
#else
const int REPETITIONS = 10;
const int N = 100000;
const int RANGE_MIN=10000;
#endif
const int VALID_NUMBER_OF_KEYS = 100;
const double EXPECTED_SUM = (REPETITIONS + 1) * N;

//! A minimal class that occupies N bytes. 
/** Defines default and copy constructor, and allows implicit operator&.
    Hides operator=. */
template<size_t N=tbb::internal::NFS_MaxLineSize>
class minimal: NoAssign {
private:
    int my_value;
    bool is_constructed;
    char pad[N-sizeof(int) - sizeof(bool)];
public:
    minimal() : NoAssign(), my_value(0) { ++construction_counter; is_constructed = true; }
    minimal( const minimal &m ) : NoAssign(), my_value(m.my_value) { ++construction_counter; is_constructed = true; }
    ~minimal() { ++destruction_counter; ASSERT(is_constructed, NULL); is_constructed = false; }
    void set_value( const int i ) { ASSERT(is_constructed, NULL); my_value = i; }
    int value( ) const { ASSERT(is_constructed, NULL); return my_value; }
};

//
// A helper class that simplifies writing the tests since minimal does not 
// define = or + operators.
//

template< typename T >
struct test_helper {
   static inline void init(T &e) { e = static_cast<T>(0); }  
   static inline void sum(T &e, const int addend ) { e += static_cast<T>(addend); }
   static inline void sum(T &e, const double addend ) { e += static_cast<T>(addend); }
   static inline void set(T &e, const int value ) { e = static_cast<T>(value); }
   static inline double get(const T &e ) { return static_cast<double>(e); }
};

template<size_t N>
struct test_helper<minimal<N> > {
   static inline void init(minimal<N> &sum) { sum.set_value( 0 ); }  
   static inline void sum(minimal<N> &sum, const int addend ) { sum.set_value( sum.value() + addend); }
   static inline void sum(minimal<N> &sum, const double addend ) { sum.set_value( sum.value() + static_cast<int>(addend)); }
   static inline void sum(minimal<N> &sum, const minimal<N> &addend ) { sum.set_value( sum.value() + addend.value()); }
   static inline void set(minimal<N> &v, const int value ) { v.set_value( static_cast<int>(value) ); }
   static inline double get(const minimal<N> &sum ) { return static_cast<double>(sum.value()); }
};

//! Tag class used to make certain constructors hard to invoke accidentally.
struct SecretTagType {} SecretTag;

//// functors and routines for initialization and combine

// Addition

template <typename T>
struct FunctorAddCombineRef {
    T operator()(const T& left, const T& right) const {
        return left+right;
    }
};

template <size_t N>
struct FunctorAddCombineRef<minimal<N> > {
    minimal<N> operator()(const minimal<N>& left, const minimal<N>& right) const {
        minimal<N> result;
        result.set_value( left.value() + right.value() ); 
        return result;
    }
};

//! Counts instances of FunctorFinit
static tbb::atomic<int> FinitCounter;

template <typename T, int Value>
struct FunctorFinit {
    FunctorFinit( const FunctorFinit& ) {++FinitCounter;}
    FunctorFinit( SecretTagType ) {++FinitCounter;}
    ~FunctorFinit() {--FinitCounter;}
    T operator()() { return Value; }
};

template <size_t N, int Value>
struct FunctorFinit<minimal<N>,Value> {
    FunctorFinit( const FunctorFinit& ) {++FinitCounter;}
    FunctorFinit( SecretTagType ) {++FinitCounter;}
    ~FunctorFinit() {--FinitCounter;}
    minimal<N> operator()() {   
        minimal<N> result;
        result.set_value( Value );
        return result;
    }
};

template <typename T>
struct FunctorAddCombine {
    T operator()(T left, T right ) const {
        return FunctorAddCombineRef<T>()( left, right );
    }
};

template <typename T>
T my_combine_ref( const T &left, const T &right) { 
    return FunctorAddCombineRef<T>()( left, right );
}

template <typename T>
T my_combine( T left, T right) { return my_combine_ref(left,right); }

template <typename T>
class combine_one_helper {
public:
    combine_one_helper(T& _result) : my_result(_result) {}
    void operator()(const T& new_bit) { test_helper<T>::sum(my_result, new_bit); }
    combine_one_helper& operator=(const combine_one_helper& other) { 
        test_helper<T>::set(my_result, test_helper<T>::get(other)); 
        return *this; 
    }
private:
    T& my_result;
};

//// end functors and routines

template< typename T >
void run_serial_scalar_tests(const char *test_name) {
    tbb::tick_count t0;
    T sum;
    test_helper<T>::init(sum);

    REMARK("Testing serial %s... ", test_name);  
    for (int t = -1; t < REPETITIONS; ++t) {
        if (Verbose && t == 0) t0 = tbb::tick_count::now(); 
        for (int i = 0; i < N; ++i) {
            test_helper<T>::sum(sum,1); 
        }
    }
 
    double result_value = test_helper<T>::get(sum);
    ASSERT( EXPECTED_SUM == result_value, NULL);
    REMARK("done\nserial %s, 0, %g, %g\n", test_name, result_value, ( tbb::tick_count::now() - t0).seconds());
}


template <typename T>
class parallel_scalar_body: NoAssign {
    
    tbb::enumerable_thread_specific<T> &sums;
 
public:

    parallel_scalar_body ( tbb::enumerable_thread_specific<T> &_sums ) : sums(_sums) { }

    void operator()( const tbb::blocked_range<int> &r ) const {
        for (int i = r.begin(); i != r.end(); ++i) 
            test_helper<T>::sum( sums.local(), 1 );
    }
   
};

template< typename T >
void run_parallel_scalar_tests_nocombine(const char *test_name) {

    typedef tbb::enumerable_thread_specific<T> ets_type;

    // We assume that static_sums zero-initialized or has a default constructor that zeros it.
    static ets_type static_sums = ets_type( T() );

    T exemplar;
    test_helper<T>::init(exemplar);
    T exemplar23;
    test_helper<T>::set(exemplar23,23);

    for (int p = MinThread; p <= MaxThread; ++p) { 
        REMARK("Testing parallel %s on %d thread(s)... ", test_name, p); 
        tbb::task_scheduler_init init(p);
        tbb::tick_count t0;

        T iterator_sum;
        test_helper<T>::init(iterator_sum);

        T finit_ets_sum;
        test_helper<T>::init(finit_ets_sum);

        T const_iterator_sum; 
        test_helper<T>::init(const_iterator_sum);

        T range_sum;
        test_helper<T>::init(range_sum);

        T const_range_sum;
        test_helper<T>::init(const_range_sum);

        T cconst_sum;
        test_helper<T>::init(cconst_sum);

        T assign_sum;
        test_helper<T>::init(assign_sum);

        T cassgn_sum;
        test_helper<T>::init(cassgn_sum);
        T non_cassgn_sum;
        test_helper<T>::init(non_cassgn_sum);

        T static_sum;
        test_helper<T>::init(static_sum);

        for (int t = -1; t < REPETITIONS; ++t) {
            if (Verbose && t == 0) t0 = tbb::tick_count::now(); 

            static_sums.clear();

            ets_type sums(exemplar);
            FunctorFinit<T,0> my_finit(SecretTag);
            ets_type finit_ets(my_finit);

            ASSERT( sums.empty(), NULL);
            tbb::parallel_for( tbb::blocked_range<int>( 0, N, RANGE_MIN ), parallel_scalar_body<T>( sums ) );
            ASSERT( !sums.empty(), NULL);

            ASSERT( finit_ets.empty(), NULL);
            tbb::parallel_for( tbb::blocked_range<int>( 0, N, RANGE_MIN ), parallel_scalar_body<T>( finit_ets ) );
            ASSERT( !finit_ets.empty(), NULL);

            ASSERT(static_sums.empty(), NULL);
            tbb::parallel_for( tbb::blocked_range<int>( 0, N, RANGE_MIN ), parallel_scalar_body<T>( static_sums ) );
            ASSERT( !static_sums.empty(), NULL);

            // use iterator
            typename ets_type::size_type size = 0;
            for ( typename ets_type::iterator i = sums.begin(); i != sums.end(); ++i ) {
                 ++size;
                 test_helper<T>::sum(iterator_sum, *i);
            }
            ASSERT( sums.size() == size, NULL);

            // use const_iterator
            for ( typename ets_type::const_iterator i = sums.begin(); i != sums.end(); ++i ) {
                 test_helper<T>::sum(const_iterator_sum, *i);
            }
           
            // use range_type
            typename ets_type::range_type r = sums.range();  
            for ( typename ets_type::range_type::const_iterator i = r.begin(); i != r.end(); ++i ) {
                 test_helper<T>::sum(range_sum, *i);
            }
           
            // use const_range_type
            typename ets_type::const_range_type cr = sums.range();  
            for ( typename ets_type::const_range_type::iterator i = cr.begin(); i != cr.end(); ++i ) {
                 test_helper<T>::sum(const_range_sum, *i);
            }

            // test copy constructor, with TLS-cached locals
            typedef typename tbb::enumerable_thread_specific<T, tbb::cache_aligned_allocator<T>, tbb::ets_key_per_instance> cached_ets_type;

            cached_ets_type cconst(sums); 

            for ( typename cached_ets_type::const_iterator i = cconst.begin(); i != cconst.end(); ++i ) {
                 test_helper<T>::sum(cconst_sum, *i);
            }
           
            // test assignment
            ets_type assigned;
            assigned = sums;

            for ( typename ets_type::const_iterator i = assigned.begin(); i != assigned.end(); ++i ) {
                 test_helper<T>::sum(assign_sum, *i);
            }

            // test assign to and from cached locals
            cached_ets_type cassgn;
            cassgn = sums;
            for ( typename cached_ets_type::const_iterator i = cassgn.begin(); i != cassgn.end(); ++i ) {
                 test_helper<T>::sum(cassgn_sum, *i);
            }

            ets_type non_cassgn;
            non_cassgn = cassgn;
            for ( typename ets_type::const_iterator i = non_cassgn.begin(); i != non_cassgn.end(); ++i ) {
                 test_helper<T>::sum(non_cassgn_sum, *i);
            }

            // test finit-initialized ets
            for(typename ets_type::const_iterator i = finit_ets.begin(); i != finit_ets.end(); ++i) {
                test_helper<T>::sum(finit_ets_sum, *i);
            }

            // test static ets
            for(typename ets_type::const_iterator i = static_sums.begin(); i != static_sums.end(); ++i) {
                test_helper<T>::sum(static_sum, *i);
            }

        }

        ASSERT( EXPECTED_SUM == test_helper<T>::get(iterator_sum), NULL);
        ASSERT( EXPECTED_SUM == test_helper<T>::get(const_iterator_sum), NULL);
        ASSERT( EXPECTED_SUM == test_helper<T>::get(range_sum), NULL);
        ASSERT( EXPECTED_SUM == test_helper<T>::get(const_range_sum), NULL);

        ASSERT( EXPECTED_SUM == test_helper<T>::get(cconst_sum), NULL);
        ASSERT( EXPECTED_SUM == test_helper<T>::get(assign_sum), NULL);
        ASSERT( EXPECTED_SUM == test_helper<T>::get(cassgn_sum), NULL);
        ASSERT( EXPECTED_SUM == test_helper<T>::get(non_cassgn_sum), NULL);
        ASSERT( EXPECTED_SUM == test_helper<T>::get(finit_ets_sum), NULL);
        ASSERT( EXPECTED_SUM == test_helper<T>::get(static_sum), NULL);

        REMARK("done\nparallel %s, %d, %g, %g\n", test_name, p, test_helper<T>::get(iterator_sum), 
                                                      ( tbb::tick_count::now() - t0).seconds());
    }
}

template< typename T >
void run_parallel_scalar_tests(const char *test_name) {

    typedef tbb::enumerable_thread_specific<T> ets_type;

    // We assume that static_sums zero-initialized or has a default constructor that zeros it.
    static ets_type static_sums = ets_type( T() );

    T exemplar;
    test_helper<T>::init(exemplar);

    run_parallel_scalar_tests_nocombine<T>(test_name);

    for (int p = MinThread; p <= MaxThread; ++p) { 
        REMARK("Testing parallel %s on %d thread(s)... ", test_name, p); 
        tbb::task_scheduler_init init(p);
        tbb::tick_count t0;

        T combine_sum;
        test_helper<T>::init(combine_sum);

        T combine_ref_sum;
        test_helper<T>::init(combine_ref_sum);

        T combine_one_sum;
        test_helper<T>::init(combine_one_sum);

        T static_sum;
        test_helper<T>::init(static_sum);

        for (int t = -1; t < REPETITIONS; ++t) {
            if (Verbose && t == 0) t0 = tbb::tick_count::now(); 

            static_sums.clear();

            ets_type sums(exemplar);

            ASSERT( sums.empty(), NULL);
            tbb::parallel_for( tbb::blocked_range<int>( 0, N, RANGE_MIN ), parallel_scalar_body<T>( sums ) );
            ASSERT( !sums.empty(), NULL);

            ASSERT(static_sums.empty(), NULL);
            tbb::parallel_for( tbb::blocked_range<int>( 0, N, RANGE_MIN ), parallel_scalar_body<T>( static_sums ) );
            ASSERT( !static_sums.empty(), NULL);


            // Use combine
            test_helper<T>::sum(combine_sum, sums.combine(my_combine<T>));
            test_helper<T>::sum(combine_ref_sum, sums.combine(my_combine_ref<T>));
            test_helper<T>::sum(static_sum, static_sums.combine(my_combine<T>));

            combine_one_helper<T> my_helper(combine_one_sum);
            sums.combine_each(my_helper);
            }


        ASSERT( EXPECTED_SUM == test_helper<T>::get(combine_sum), NULL);
        ASSERT( EXPECTED_SUM == test_helper<T>::get(combine_ref_sum), NULL);
        ASSERT( EXPECTED_SUM == test_helper<T>::get(static_sum), NULL);

        REMARK("done\nparallel combine %s, %d, %g, %g\n", test_name, p, test_helper<T>::get(combine_sum), 
                                                      ( tbb::tick_count::now() - t0).seconds());
    }
}

template <typename T>
class parallel_vector_for_body: NoAssign {
    
    tbb::enumerable_thread_specific< std::vector<T, tbb::tbb_allocator<T> > > &locals;
 
public:

    parallel_vector_for_body ( tbb::enumerable_thread_specific< std::vector<T, tbb::tbb_allocator<T> > > &_locals ) : locals(_locals) { }

    void operator()( const tbb::blocked_range<int> &r ) const {
        T one;
        test_helper<T>::set(one, 1);

        for (int i = r.begin(); i < r.end(); ++i) {
            locals.local().push_back( one );
        }
    }
   
};

template <typename R, typename T>
struct parallel_vector_reduce_body {

    T sum;    
    size_t count;    

    parallel_vector_reduce_body ( ) : count(0) { test_helper<T>::init(sum); }
    parallel_vector_reduce_body ( parallel_vector_reduce_body<R, T> &, tbb::split ) : count(0) {  test_helper<T>::init(sum); }

    void operator()( const R &r ) {
        for (typename R::iterator ri = r.begin(); ri != r.end(); ++ri) {
            const std::vector< T, tbb::tbb_allocator<T>  > &v = *ri; 
            ++count;
            for (typename std::vector<T, tbb::tbb_allocator<T> >::const_iterator vi = v.begin(); vi != v.end(); ++vi) {
                test_helper<T>::sum(sum, *vi);
            }
        }
    }

    void join( const parallel_vector_reduce_body &b ) {
        test_helper<T>::sum(sum,b.sum);
        count += b.count;
    }
   
};

template< typename T >
void run_parallel_vector_tests(const char *test_name) {
    tbb::tick_count t0;
    typedef std::vector<T, tbb::tbb_allocator<T> > container_type;

    for (int p = MinThread; p <= MaxThread; ++p) {
        REMARK("Testing parallel %s on %d thread(s)... ", test_name, p);
        tbb::task_scheduler_init init(p);

        T sum;
        test_helper<T>::init(sum);

        for (int t = -1; t < REPETITIONS; ++t) {
            if (Verbose && t == 0) t0 = tbb::tick_count::now(); 
            typedef typename tbb::enumerable_thread_specific< container_type > ets_type;
            ets_type vs;

            ASSERT( vs.empty(), NULL);
            tbb::parallel_for ( tbb::blocked_range<int> (0, N, RANGE_MIN), parallel_vector_for_body<T>( vs ) );
            ASSERT( !vs.empty(), NULL);

            // copy construct
            ets_type vs2(vs); // this causes an assertion failure, related to allocators...

            // assign
            ets_type vs3;
            vs3 = vs;

            parallel_vector_reduce_body< typename tbb::enumerable_thread_specific< std::vector< T, tbb::tbb_allocator<T>  > >::const_range_type, T > pvrb;
            tbb::parallel_reduce ( vs.range(1), pvrb );

            test_helper<T>::sum(sum, pvrb.sum);

            ASSERT( vs.size() == pvrb.count, NULL);

            tbb::flattened2d<ets_type> fvs = flatten2d(vs);
            size_t ccount = fvs.size();
            size_t elem_cnt = 0;
            for(typename tbb::flattened2d<ets_type>::const_iterator i = fvs.begin(); i != fvs.end(); ++i) {
                ++elem_cnt;
            };
            ASSERT(ccount == elem_cnt, NULL);

            elem_cnt = 0;
            for(typename tbb::flattened2d<ets_type>::iterator i = fvs.begin(); i != fvs.end(); ++i) {
                ++elem_cnt;
            };
            ASSERT(ccount == elem_cnt, NULL);
        }

        double result_value = test_helper<T>::get(sum);
        ASSERT( EXPECTED_SUM == result_value, NULL);
        REMARK("done\nparallel %s, %d, %g, %g\n", test_name, p, result_value, ( tbb::tick_count::now() - t0).seconds());
    }
}

template<typename T>
void run_cross_type_vector_tests(const char *test_name) {
    tbb::tick_count t0;
    typedef std::vector<T, tbb::tbb_allocator<T> > container_type;

    for (int p = MinThread; p <= MaxThread; ++p) {
        REMARK("Testing parallel %s on %d thread(s)... ", test_name, p);
        tbb::task_scheduler_init init(p);

        T sum;
        test_helper<T>::init(sum);

        for (int t = -1; t < REPETITIONS; ++t) {
            if (Verbose && t == 0) t0 = tbb::tick_count::now(); 
            typedef typename tbb::enumerable_thread_specific< container_type, tbb::cache_aligned_allocator<container_type>, tbb::ets_no_key > ets_nokey_type;
            typedef typename tbb::enumerable_thread_specific< container_type, tbb::cache_aligned_allocator<container_type>, tbb::ets_key_per_instance > ets_tlskey_type;
            ets_nokey_type vs;

            ASSERT( vs.empty(), NULL);
            tbb::parallel_for ( tbb::blocked_range<int> (0, N, RANGE_MIN), parallel_vector_for_body<T>( vs ) );
            ASSERT( !vs.empty(), NULL);

            // copy construct
            ets_tlskey_type vs2(vs);

            // assign
            ets_nokey_type vs3;
            vs3 = vs2;

            parallel_vector_reduce_body< typename tbb::enumerable_thread_specific< std::vector< T, tbb::tbb_allocator<T>  > >::const_range_type, T > pvrb;
            tbb::parallel_reduce ( vs3.range(1), pvrb );

            test_helper<T>::sum(sum, pvrb.sum);

            ASSERT( vs3.size() == pvrb.count, NULL);

            tbb::flattened2d<ets_nokey_type> fvs = flatten2d(vs3);
            size_t ccount = fvs.size();
            size_t elem_cnt = 0;
            for(typename tbb::flattened2d<ets_nokey_type>::const_iterator i = fvs.begin(); i != fvs.end(); ++i) {
                ++elem_cnt;
            };
            ASSERT(ccount == elem_cnt, NULL);

            elem_cnt = 0;
            for(typename tbb::flattened2d<ets_nokey_type>::iterator i = fvs.begin(); i != fvs.end(); ++i) {
                ++elem_cnt;
            };
            ASSERT(ccount == elem_cnt, NULL);
        }

        double result_value = test_helper<T>::get(sum);
        ASSERT( EXPECTED_SUM == result_value, NULL);
        REMARK("done\nparallel %s, %d, %g, %g\n", test_name, p, result_value, ( tbb::tick_count::now() - t0).seconds());
    }
}

template< typename T >
void run_serial_vector_tests(const char *test_name) {
    tbb::tick_count t0;
    T sum;
    test_helper<T>::init(sum);
    T one;
    test_helper<T>::set(one, 1);

    REMARK("Testing serial %s... ", test_name);
    for (int t = -1; t < REPETITIONS; ++t) {
        if (Verbose && t == 0) t0 = tbb::tick_count::now(); 
        std::vector<T, tbb::tbb_allocator<T> > v; 
        for (int i = 0; i < N; ++i) {
            v.push_back( one );
        }
        for (typename std::vector<T, tbb::tbb_allocator<T> >::const_iterator i = v.begin(); i != v.end(); ++i) 
            test_helper<T>::sum(sum, *i); 
    }

    double result_value = test_helper<T>::get(sum);
    ASSERT( EXPECTED_SUM == result_value, NULL);
    REMARK("done\nserial %s, 0, %g, %g\n", test_name, result_value, ( tbb::tick_count::now() - t0).seconds());
}

const size_t line_size = tbb::internal::NFS_MaxLineSize;

void 
run_serial_tests() {
    run_serial_scalar_tests<int>("int");
    run_serial_scalar_tests<double>("double");
    run_serial_scalar_tests<minimal<> >("minimal<>");
    run_serial_vector_tests<int>("std::vector<int, tbb::tbb_allocator<int> >");
    run_serial_vector_tests<double>("std::vector<double, tbb::tbb_allocator<double> >");
}

void 
run_parallel_tests() {
    run_parallel_scalar_tests<int>("int");
    run_parallel_scalar_tests<double>("double");
    run_parallel_scalar_tests_nocombine<minimal<> >("minimal<>");
    run_parallel_vector_tests<int>("std::vector<int, tbb::tbb_allocator<int> >");
    run_parallel_vector_tests<double>("std::vector<double, tbb::tbb_allocator<double> >");
}

void
run_cross_type_tests() {
    // cross-type scalar tests are part of run_serial_scalar_tests
    run_cross_type_vector_tests<int>("std::vector<int, tbb::tbb_allocator<int> >");
    run_parallel_vector_tests<double>("std::vector<double, tbb::tbb_allocator<double> >");
}

typedef tbb::enumerable_thread_specific<minimal<line_size> > flogged_ets;

class set_body {
    flogged_ets *a;

public:
    set_body( flogged_ets*_a ) : a(_a) { }

    void operator() ( ) const {
        for (int i = 0; i < VALID_NUMBER_OF_KEYS; ++i) {
            a[i].local().set_value(i + 1);
        }
    }
 
};

void do_tbb_threads( int max_threads, flogged_ets a[] ) {
    std::vector< tbb::tbb_thread * > threads;

    for (int p = 0; p < max_threads; ++p) { 
        threads.push_back( new tbb::tbb_thread ( set_body( a ) ) ); 
    }

    for (int p = 0; p < max_threads; ++p) {
        threads[p]->join();
    }

    for(int p = 0; p < max_threads; ++p) {
        delete threads[p];
    }
}

void
flog_key_creation_and_deletion() {
    const int FLOG_REPETITIONS = 100;

    for (int p = MinThread; p <= MaxThread; ++p) { 
        REMARK("Testing repeated deletes on %d threads... ", p);

        for (int j = 0; j < FLOG_REPETITIONS; ++j) {
            construction_counter = 0;
            destruction_counter = 0;

            // causes VALID_NUMER_OF_KEYS exemplar instances to be constructed 
            flogged_ets* a = new flogged_ets[VALID_NUMBER_OF_KEYS];
            ASSERT(int(construction_counter) == 0, NULL);   // no exemplars or actual locals have been constructed
            ASSERT(int(destruction_counter) == 0, NULL);    // and none have been destroyed

            // causes p * VALID_NUMBER_OF_KEYS minimals to be created
            do_tbb_threads(p, a); 

            for (int i = 0; i < VALID_NUMBER_OF_KEYS; ++i) {
                int pcnt = 0;
                for ( flogged_ets::iterator tli = a[i].begin(); tli != a[i].end(); ++tli ) {
                    ASSERT( (*tli).value() == i+1, NULL );
                    ++pcnt;
                }
                ASSERT( pcnt == p, NULL);  // should be one local per thread.
            }
            delete[] a;
        }

        ASSERT( int(construction_counter) == (p)*VALID_NUMBER_OF_KEYS, NULL );
        ASSERT( int(destruction_counter) == (p)*VALID_NUMBER_OF_KEYS, NULL );

        REMARK("done\nTesting repeated clears on %d threads... ", p);

        construction_counter = 0;
        destruction_counter = 0;

        // causes VALID_NUMER_OF_KEYS exemplar instances to be constructed 
        flogged_ets* a = new flogged_ets[VALID_NUMBER_OF_KEYS];
 
        for (int j = 0; j < FLOG_REPETITIONS; ++j) {

            // causes p * VALID_NUMBER_OF_KEYS minimals to be created
            do_tbb_threads(p, a);

            for (int i = 0; i < VALID_NUMBER_OF_KEYS; ++i) {
                for ( flogged_ets::iterator tli = a[i].begin(); tli != a[i].end(); ++tli ) {
                    ASSERT( (*tli).value() == i+1, NULL );
                }
                a[i].clear();
                ASSERT( static_cast<int>(a[i].end() - a[i].begin()) == 0, NULL );
            }

        }

        delete[] a;

        ASSERT( int(construction_counter) == (FLOG_REPETITIONS*p)*VALID_NUMBER_OF_KEYS, NULL );
        ASSERT( int(destruction_counter) == (FLOG_REPETITIONS*p)*VALID_NUMBER_OF_KEYS, NULL );

        REMARK("done\n");
    }

}

template <typename inner_container>
void 
flog_segmented_interator() {

    bool found_error = false;
    typedef typename inner_container::value_type T;
    typedef std::vector< inner_container > nested_vec;
    inner_container my_inner_container;
    my_inner_container.clear();
    nested_vec my_vec;

    // simple nested vector (neither level empty)
    const int maxval = 10;
    for(int i=0; i < maxval; i++) {
        my_vec.push_back(my_inner_container);
        for(int j = 0; j < maxval; j++) {
            my_vec.at(i).push_back((T)(maxval * i + j));
        }
    }

    tbb::internal::segmented_iterator<nested_vec, T> my_si(my_vec);

    T ii;
    for(my_si=my_vec.begin(), ii=0; my_si != my_vec.end(); ++my_si, ++ii) {
        if((*my_si) != ii) {
            found_error = true;
            REMARK( "*my_si=%d\n", int(*my_si));
        }
    }

    // outer level empty
    my_vec.clear();
    for(my_si=my_vec.begin(); my_si != my_vec.end(); ++my_si) {
        found_error = true;
    }

    // inner levels empty
    my_vec.clear();
    for(int i =0; i < maxval; ++i) {
        my_vec.push_back(my_inner_container);
    }
    for(my_si = my_vec.begin(); my_si != my_vec.end(); ++my_si) {
        found_error = true;
    }

    // every other inner container is empty
    my_vec.clear();
    for(int i=0; i < maxval; ++i) {
        my_vec.push_back(my_inner_container);
        if(i%2) {
            for(int j = 0; j < maxval; ++j) {
                my_vec.at(i).push_back((T)(maxval * (i/2) + j));
            }
        }
    }
    for(my_si = my_vec.begin(), ii=0; my_si != my_vec.end(); ++my_si, ++ii) {
        if((*my_si) != ii) {
            found_error = true;
            REMARK("*my_si=%d, ii=%d\n", (int)(*my_si), (int)ii);
        }
    }

    tbb::internal::segmented_iterator<nested_vec, const T> my_csi(my_vec);
    for(my_csi=my_vec.begin(), ii=0; my_csi != my_vec.end(); ++my_csi, ++ii) {
        if((*my_csi) != ii) {
            found_error = true;
            REMARK( "*my_csi=%d\n", int(*my_csi));
        }
    }

    // outer level empty
    my_vec.clear();
    for(my_csi=my_vec.begin(); my_csi != my_vec.end(); ++my_csi) {
        found_error = true;
    }

    // inner levels empty
    my_vec.clear();
    for(int i =0; i < maxval; ++i) {
        my_vec.push_back(my_inner_container);
    }
    for(my_csi = my_vec.begin(); my_csi != my_vec.end(); ++my_csi) {
        found_error = true;
    }

    // every other inner container is empty
    my_vec.clear();
    for(int i=0; i < maxval; ++i) {
        my_vec.push_back(my_inner_container);
        if(i%2) {
            for(int j = 0; j < maxval; ++j) {
                my_vec.at(i).push_back((T)(maxval * (i/2) + j));
            }
        }
    }
    for(my_csi = my_vec.begin(), ii=0; my_csi != my_vec.end(); ++my_csi, ++ii) {
        if((*my_csi) != ii) {
            found_error = true;
            REMARK("*my_csi=%d, ii=%d\n", (int)(*my_csi), (int)ii);
        }
    }


    if(found_error) REPORT("segmented_iterator failed\n");
}

template <typename Key, typename Val>
void
flog_segmented_iterator_map() {
   typedef typename std::map<Key, Val> my_map;
   typedef std::vector< my_map > nested_vec;
   my_map my_inner_container;
   my_inner_container.clear();
   nested_vec my_vec;
   my_vec.clear();
   bool found_error = false;

   // simple nested vector (neither level empty)
   const int maxval = 4;
   for(int i=0; i < maxval; i++) {
       my_vec.push_back(my_inner_container);
       for(int j = 0; j < maxval; j++) {
           my_vec.at(i).insert(std::make_pair<Key,Val>(maxval * i + j, 2*(maxval*i + j)));
       }
   }

   tbb::internal::segmented_iterator<nested_vec, std::pair<const Key, Val> > my_si(my_vec);
   Key ii;
   for(my_si=my_vec.begin(), ii=0; my_si != my_vec.end(); ++my_si, ++ii) {
       if(((*my_si).first != ii) || ((*my_si).second != 2*ii)) {
           found_error = true;
           REMARK( "ii=%d, (*my_si).first=%d, second=%d\n",ii, int((*my_si).first), int((*my_si).second));
       }
   }

   tbb::internal::segmented_iterator<nested_vec, const std::pair<const Key, Val> > my_csi(my_vec);
   for(my_csi=my_vec.begin(), ii=0; my_csi != my_vec.end(); ++my_csi, ++ii) {
       if(((*my_csi).first != ii) || ((*my_csi).second != 2*ii)) {
           found_error = true;
           REMARK( "ii=%d, (*my_csi).first=%d, second=%d\n",ii, int((*my_csi).first), int((*my_csi).second));
       }
   }
   if(found_error) REPORT("segmented_iterator_map failed\n");
}

void
run_segmented_iterator_tests() {
   // only the following containers can be used with the segmented iterator.
   REMARK("Running Segmented Iterator Tests\n");
   flog_segmented_interator<std::vector< int > >();
   flog_segmented_interator<std::vector< double > >();
   flog_segmented_interator<std::deque< int > >();
   flog_segmented_interator<std::deque< double > >();
   flog_segmented_interator<std::list< int > >();
   flog_segmented_interator<std::list< double > >();

   flog_segmented_iterator_map<int, int>();
   flog_segmented_iterator_map<int, double>(); 
}

template <typename T>
void
run_assign_and_copy_constructor_test(const char *test_name) {
    REMARK("Testing assignment and copy construction for %s\n", test_name);

    // test initializer with exemplar
    T initializer0;
    test_helper<T>::init(initializer0);
    T initializer7;
    test_helper<T>::set(initializer7,7);
    tbb::enumerable_thread_specific<T> create1(initializer7);
    (void) create1.local();  // create an initialized value
    ASSERT(7 == test_helper<T>::get(create1.local()), NULL);

    // test copy construction with exemplar initializer
    create1.clear();
    tbb::enumerable_thread_specific<T> copy1(create1);
    (void) copy1.local();
    ASSERT(7 == test_helper<T>::get(copy1.local()), NULL);

    // test copy assignment with exemplar initializer
    create1.clear();
    tbb::enumerable_thread_specific<T> assign1(initializer0);
    assign1 = create1;
    (void) assign1.local();
    ASSERT(7 == test_helper<T>::get(assign1.local()), NULL);

    // test creation with finit function
    FunctorFinit<T,7> my_finit7(SecretTag);
    tbb::enumerable_thread_specific<T> create2(my_finit7);
    (void) create2.local();
    ASSERT(7 == test_helper<T>::get(create2.local()), NULL);

    // test copy construction with function initializer
    create2.clear();
    tbb::enumerable_thread_specific<T> copy2(create2);
    (void) copy2.local();
    ASSERT(7 == test_helper<T>::get(copy2.local()), NULL);

    // test copy assignment with function initializer
    create2.clear();
    FunctorFinit<T,0> my_finit(SecretTag);
    tbb::enumerable_thread_specific<T> assign2(my_finit);
    assign2 = create2;
    (void) assign2.local();
    ASSERT(7 == test_helper<T>::get(assign2.local()), NULL);
}

void
run_assignment_and_copy_constructor_tests() {
    REMARK("Running assignment and copy constructor tests\n");
    run_assign_and_copy_constructor_test<int>("int");
    run_assign_and_copy_constructor_test<double>("double");
    // Try class sizes that are close to a cache line in size, in order to check padding calculations.
    run_assign_and_copy_constructor_test<minimal<line_size-1> >("minimal<line_size-1>");
    run_assign_and_copy_constructor_test<minimal<line_size> >("minimal<line_size>");
    run_assign_and_copy_constructor_test<minimal<line_size+1> >("minimal<line_size+1>");
    ASSERT(FinitCounter==0, NULL);
}

// Class with no default constructor
class HasNoDefaultConstructor {
    HasNoDefaultConstructor();
public:
    HasNoDefaultConstructor( SecretTagType ) {}
};

// Initialization functor for a HasNoDefaultConstructor
struct HasNoDefaultConstructorFinit {
    HasNoDefaultConstructor operator()() {
        return HasNoDefaultConstructor(SecretTag);
    }
};

struct HasNoDefaultConstructorCombine {
    HasNoDefaultConstructor operator()( HasNoDefaultConstructor, HasNoDefaultConstructor ) {
        return HasNoDefaultConstructor(SecretTag);
    }
};

//! Test situations where only default constructor or copy constructor is required.
void TestInstantiation() {
    // Test instantiation is possible when copy constructor is not required.
    tbb::enumerable_thread_specific<NoCopy> ets1;

    // Test instantiation when default constructor is not required, because exemplar is provided.
    HasNoDefaultConstructor x(SecretTag);
    tbb::enumerable_thread_specific<HasNoDefaultConstructor> ets2(x);
    ets2.combine(HasNoDefaultConstructorCombine());

    // Test instantiation when default constructor is not required, because init function is provided.
    HasNoDefaultConstructorFinit f;
    tbb::enumerable_thread_specific<HasNoDefaultConstructor> ets3(f);
    ets3.combine(HasNoDefaultConstructorCombine());
}

class BigType {
public:
    BigType() { /* avoid cl warning C4345 about default initialization of POD types */ }
    char my_data[12 * 1024 * 1024];
};

void TestConstructor() {
    typedef tbb::enumerable_thread_specific<BigType> CounterBigType;
    // Test default constructor
    CounterBigType MyCounters; 
    // Create a local instance.
    CounterBigType::reference my_local = MyCounters.local();
    my_local.my_data[0] = 'a';
    // Test copy constructor
    CounterBigType MyCounters2(MyCounters);
}

int TestMain () {
    TestInstantiation();
    run_segmented_iterator_tests();
    flog_key_creation_and_deletion();

    if (MinThread == 0) {
        run_serial_tests();
        MinThread = 1;
    }
    if (MaxThread > 0) {
        run_parallel_tests();
        run_cross_type_tests();
    }

    run_assignment_and_copy_constructor_tests();

    TestConstructor();

    return Harness::Done;
}
