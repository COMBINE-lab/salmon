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

#define __TBB_EXTRA_DEBUG 1
#include "tbb/concurrent_unordered_map.h"
#if __TBB_INITIALIZER_LISTS_PRESENT
// These operator== are used implicitly in  test_initializer_list.h.
// For some unknown reason clang is not able to find the if they a declared after the
// inclusion of test_initializer_list.h.
template<typename container_type>
bool equal_containers( container_type const& lhs, container_type const& rhs );
template<typename Key, typename Value>
bool operator==( tbb::concurrent_unordered_map<Key, Value> const& lhs, tbb::concurrent_unordered_map<Key, Value> const& rhs ) {
    return equal_containers( lhs, rhs );
}
template<typename Key, typename Value>
bool operator==( tbb::concurrent_unordered_multimap<Key, Value> const& lhs, tbb::concurrent_unordered_multimap<Key, Value> const& rhs ) {
    return equal_containers( lhs, rhs );
}
#endif /* __TBB_INITIALIZER_LISTS_PRESENT */
#include "test_concurrent_unordered_common.h"

typedef tbb::concurrent_unordered_map<int, int, tbb::tbb_hash<int>, std::equal_to<int>, MyAllocator> MyMap;
typedef tbb::concurrent_unordered_map<int, check_type<int>, tbb::tbb_hash<int>, std::equal_to<int>, MyAllocator> MyCheckedMap;
typedef tbb::concurrent_unordered_multimap<int, int, tbb::tbb_hash<int>, std::equal_to<int>, MyAllocator> MyMultiMap;
typedef tbb::concurrent_unordered_multimap<int, check_type<int>, tbb::tbb_hash<int>, std::equal_to<int>, MyAllocator> MyCheckedMultiMap;

template <>
struct SpecialTests <MyMap> {
    static void Test( const char *str ) {
        MyMap cont( 0 );
        const MyMap &ccont( cont );

        // mapped_type& operator[](const key_type& k);
        cont[1] = 2;

        // bool empty() const;    
        ASSERT( !ccont.empty( ), "Concurrent container empty after adding an element" );

        // size_type size() const;
        ASSERT( ccont.size( ) == 1, "Concurrent container size incorrect" );

        ASSERT( cont[1] == 2, "Concurrent container value incorrect" );

        // mapped_type& at( const key_type& k );
        // const mapped_type& at(const key_type& k) const;
        ASSERT( cont.at( 1 ) == 2, "Concurrent container value incorrect" );
        ASSERT( ccont.at( 1 ) == 2, "Concurrent container value incorrect" );

        // iterator find(const key_type& k);
        MyMap::const_iterator it = cont.find( 1 );
        ASSERT( it != cont.end( ) && Value<MyMap>::get( *(it) ) == 2, "Element with key 1 not properly found" );
        cont.unsafe_erase( it );
        it = cont.find( 1 );
        ASSERT( it == cont.end( ), "Element with key 1 not properly erased" );

        REMARK( "passed -- specialized %s tests\n", str );
    }
};

template<>
class AssignBody<MyCheckedMap>: NoAssign{
    MyCheckedMap &table;
public:
    AssignBody( MyCheckedMap &t ) : NoAssign( ), table( t ) {}
    void operator()( int i ) const {
        table.insert( MyCheckedMap::value_type( i, check_type<int>( i ) ) );
    }
};

// for multimap insert (i%3)+1 items  [i,3*i], [i,3*i+1] ..
template<>
class AssignBody<MyMultiMap>: NoAssign{
    MyMultiMap &table;
public:
    AssignBody( MyMultiMap &t ) : NoAssign( ), table( t ) {}
    void operator()( int i ) const {
        for ( int j = 0; j < (i % 3) + 1; ++j ) {
            table.insert( std::pair<int, int>( i, 3 * i + j - 1 ) );
        }
    }
};

// for multimap insert (i%3)+1 items  [i,3*i], [i,3*i+1] ..
template<>
class AssignBody<MyCheckedMultiMap>: NoAssign{
    MyCheckedMultiMap &table;
public:
    AssignBody( MyCheckedMultiMap &t ) : NoAssign( ), table( t ) {}
    void operator()( int i ) const {
        for ( int j = 0; j < (i % 3) + 1; ++j ) {
            table.insert( std::pair<int, int>( i, 3 * i + j - 1 ) );
        }
    }
};

// test assumes the unordered multimap puts items in ascending order, and the insertions
// occur at the end of a range. This assumption may not always be valid.
template <>
struct SpecialTests <MyMultiMap> {
#define VALUE1 7
#define VALUE2 2
    static void Test( const char *str ) {
        MyMultiMap cont( 0 );
        const MyMultiMap &ccont( cont );
        // mapped_type& operator[](const key_type& k);
        cont.insert( std::make_pair( 1, VALUE1 ) );

        // bool empty() const;    
        ASSERT( !ccont.empty( ), "Concurrent container empty after adding an element" );

        // size_type size() const;
        ASSERT( ccont.size( ) == 1, "Concurrent container size incorrect" );
        ASSERT( (*(cont.begin( ))).second == VALUE1, "Concurrent container value incorrect" );
        ASSERT( (*(cont.equal_range( 1 )).first).second == VALUE1, "Improper value from equal_range" );
        ASSERT( (cont.equal_range( 1 )).second == cont.end( ), "Improper iterator from equal_range" );

        cont.insert( std::make_pair( 1, VALUE2 ) );

        // bool empty() const;    
        ASSERT( !ccont.empty( ), "Concurrent container empty after adding an element" );

        // size_type size() const;
        ASSERT( ccont.size( ) == 2, "Concurrent container size incorrect" );
        ASSERT( (*(cont.begin( ))).second == VALUE1, "Concurrent container value incorrect" );
        ASSERT( (*(cont.equal_range( 1 )).first).second == VALUE1, "Improper value from equal_range" );
        ASSERT( (cont.equal_range( 1 )).second == cont.end( ), "Improper iterator from equal_range" );

        // check that the second value is part of the range.
        // though I am not sure there are guarantees what order the insertions appear in the range
        // if the order differs the ASSERT above will fail already.
        std::pair<MyMultiMap::iterator, MyMultiMap::iterator> range = cont.equal_range( 1 );
        MyMultiMap::iterator ii = range.first;
        ++ii;
        ASSERT( (*ii).second == VALUE2, "Improper value for second insertion" );

        cont.insert( std::make_pair( 0, 4 ) );

        // bool empty() const;    
        ASSERT( !ccont.empty( ), "Concurrent container empty after adding an element" );

        // size_type size() const;
        ASSERT( ccont.size( ) == 3, "Concurrent container size incorrect" );
        ASSERT( (*(cont.begin( ))).second == 4, "Concurrent container value incorrect" );
        ASSERT( (*(cont.equal_range( 1 )).first).second == VALUE1, "Improper value from equal_range" );
        ASSERT( (cont.equal_range( 1 )).second == cont.end( ), "Improper iterator from equal_range" );

        REMARK( "passed -- specialized %s tests\n", str );
    }
};

#if __TBB_RANGE_BASED_FOR_PRESENT
#include "test_range_based_for.h"
// Add the similar test for concurrent_unordered_set.
void TestRangeBasedFor() {
    using namespace range_based_for_support_tests;

    REMARK( "testing range based for loop compatibility \n" );
    typedef tbb::concurrent_unordered_map<int, int> cu_map;
    cu_map a_cu_map;
    const int sequence_length = 100;
    for ( int i = 1; i <= sequence_length; ++i ) {
        a_cu_map.insert( cu_map::value_type( i, i ) );
    }

    ASSERT( range_based_for_accumulate( a_cu_map, pair_second_summer(), 0 ) == gauss_summ_of_int_sequence( sequence_length ), "incorrect accumulated value generated via range based for ?" );
}
#endif //if __TBB_RANGE_BASED_FOR_PRESENT

#if __TBB_CPP11_RVALUE_REF_PRESENT
struct cu_map_type : unordered_move_traits_base {
    template<typename element_type, typename allocator_type>
    struct apply {
        typedef tbb::concurrent_unordered_map<element_type, element_type, tbb::tbb_hash<element_type>, std::equal_to<element_type>, allocator_type > type;
    };

    typedef FooPairIterator init_iterator_type;
};

struct cu_multimap_type : unordered_move_traits_base {
    template<typename element_type, typename allocator_type>
    struct apply {
        typedef tbb::concurrent_unordered_multimap<element_type, element_type, tbb::tbb_hash<element_type>, std::equal_to<element_type>, allocator_type > type;
    };

    typedef FooPairIterator init_iterator_type;
};
#endif /* __TBB_CPP11_RVALUE_REF_PRESENT */

template <typename Table>
class TestOperatorSquareBrackets : NoAssign {
    typedef typename Table::value_type ValueType;
    Table &my_c;
    const ValueType &my_value;
public:
    TestOperatorSquareBrackets( Table &c, const ValueType &value ) : my_c( c ), my_value( value ) {}
    void operator()() const {
        ASSERT( Harness::IsEqual()(my_c[my_value.first], my_value.second), NULL );
    }
};

template <bool defCtorPresent, typename Key, typename Element, typename Hasher, typename Equality, typename Allocator>
void TestMapSpecificMethods( tbb::concurrent_unordered_map<Key, Element, Hasher, Equality, Allocator> &c,
    const typename tbb::concurrent_unordered_map<Key, Element, Hasher, Equality, Allocator>::value_type &value ) {
    typedef tbb::concurrent_unordered_map<Key, Element, Hasher, Equality, Allocator> Table;
    CallIf<defCtorPresent>()(TestOperatorSquareBrackets<Table>( c, value ));
    ASSERT( Harness::IsEqual()(c.at( value.first ), value.second), NULL );
    const Table &constC = c;
    ASSERT( Harness::IsEqual()(constC.at( value.first ), value.second), NULL );
}

template <bool defCtorPresent, typename ValueType>
void TestTypesMap( const std::list<ValueType> &lst ) {
    typedef typename ValueType::first_type KeyType;
    typedef typename ValueType::second_type ElemType;
    TypeTester< defCtorPresent, tbb::concurrent_unordered_map<KeyType, ElemType, tbb::tbb_hash<KeyType>, Harness::IsEqual>,
        tbb::concurrent_unordered_map< KeyType, ElemType, tbb::tbb_hash<KeyType>, Harness::IsEqual, debug_allocator<ValueType> > >( lst );
    TypeTester< defCtorPresent, tbb::concurrent_unordered_multimap<KeyType, ElemType, tbb::tbb_hash<KeyType>, Harness::IsEqual>,
        tbb::concurrent_unordered_multimap< KeyType, ElemType, tbb::tbb_hash<KeyType>, Harness::IsEqual, debug_allocator<ValueType> > >( lst );
}

void TestTypes() {
    const int NUMBER = 10;

    std::list< std::pair<const int, int> > arrIntInt;
    for ( int i = 0; i < NUMBER; ++i ) arrIntInt.push_back( std::make_pair( i, NUMBER - i ) );
    TestTypesMap</*def_ctor_present = */true>( arrIntInt );

    std::list< std::pair< const int, tbb::atomic<int> > > arrIntTbb;
    for ( int i = 0; i < NUMBER; ++i ) {
        tbb::atomic<int> b;
        b = NUMBER - i;
        arrIntTbb.push_back( std::make_pair( i, b ) );
    }
    TestTypesMap</*defCtorPresent = */true>( arrIntTbb );

#if __TBB_CPP11_REFERENCE_WRAPPER_PRESENT
    std::list< std::pair<const std::reference_wrapper<const int>, int> > arrRefInt;
    for ( std::list< std::pair<const int, int> >::iterator it = arrIntInt.begin(); it != arrIntInt.end(); ++it )
        arrRefInt.push_back( std::make_pair( std::reference_wrapper<const int>( it->first ), it->second ) );
    TestTypesMap</*defCtorPresent = */true>( arrRefInt );

    std::list< std::pair<const int, std::reference_wrapper<int> > > arrIntRef;
    for ( std::list< std::pair<const int, int> >::iterator it = arrIntInt.begin(); it != arrIntInt.end(); ++it ) {
        // Using std::make_pair below causes compilation issues with early implementations of std::reference_wrapper.
        arrIntRef.push_back( std::pair<const int, std::reference_wrapper<int> >( it->first, std::reference_wrapper<int>( it->second ) ) );
    }
    TestTypesMap</*defCtorPresent = */false>( arrIntRef );
#endif /* __TBB_CPP11_REFERENCE_WRAPPER_PRESENT */

#if __TBB_CPP11_SMART_POINTERS_PRESENT
    std::list< std::pair< const std::shared_ptr<int>, std::shared_ptr<int> > > arrShrShr;
    for ( int i = 0; i < NUMBER; ++i ) arrShrShr.push_back( std::make_pair( std::make_shared<int>( i ), std::make_shared<int>( NUMBER - i ) ) );
    TestTypesMap</*defCtorPresent = */true>( arrShrShr );

    std::list< std::pair< const std::weak_ptr<int>, std::weak_ptr<int> > > arrWkWk;
    std::copy( arrShrShr.begin(), arrShrShr.end(), std::back_inserter( arrWkWk ) );
    TestTypesMap</*defCtorPresent = */true>( arrWkWk );
#endif /* __TBB_CPP11_SMART_POINTERS_PRESENT */
}

int TestMain() {
    test_machine();

    test_basic<MyMap>( "concurrent unordered Map" );
    test_concurrent<MyMap>( "concurrent unordered Map" );
    test_basic<MyMultiMap>( "concurrent unordered MultiMap" );
    test_concurrent<MyMultiMap>( "concurrent unordered MultiMap" );
    test_concurrent<MyMultiMap>( "concurrent unordered MultiMap asymptotic", true );

    { Check<MyCheckedMap::value_type> checkit; test_basic<MyCheckedMap>( "concurrent unordered map (checked)" ); }
    { Check<MyCheckedMap::value_type> checkit; test_concurrent<MyCheckedMap>( "concurrent unordered map (checked)" ); }

    { Check<MyCheckedMultiMap::value_type> checkit; test_basic<MyCheckedMultiMap>( "concurrent unordered MultiMap (checked)" ); }
    { Check<MyCheckedMultiMap::value_type> checkit; test_concurrent<MyCheckedMultiMap>( "concurrent unordered MultiMap (checked)" ); }

#if __TBB_INITIALIZER_LISTS_PRESENT
    TestInitList< tbb::concurrent_unordered_map<int, int>,
                  tbb::concurrent_unordered_multimap<int, int> >( {{1,1},{2,2},{3,3},{4,4},{5,5}} );
#endif /* __TBB_INITIALIZER_LISTS_PRESENT */

#if __TBB_RANGE_BASED_FOR_PRESENT
    TestRangeBasedFor();
#endif

#if __TBB_CPP11_RVALUE_REF_PRESENT
    test_rvalue_ref_support<cu_map_type>( "concurrent unordered map" );
    test_rvalue_ref_support<cu_multimap_type>( "concurrent unordered multiset" );
#endif //__TBB_CPP11_RVALUE_REF_PRESENT
    
    TestTypes();

    return Harness::Done;
}
