// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

// Catch includes
#include <catch2/catch.hpp>

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Types.hpp>
using namespace Reaktoro;

TEST_CASE("Testing Algorithms", "[Algorithms]")
{
    Vec<int> nums = { 1, 2, 3, 4, 5, 6, 7, 8, 9};
    Vec<int> res;
    Strings strs;

    //-------------------------------------------------------------------------
    // TESTING METHOD: index
    //-------------------------------------------------------------------------
    for(auto n : nums)
        REQUIRE( index(nums, n) == n - 1);

    REQUIRE( index(nums, 1000) >= nums.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: indexfn
    //-------------------------------------------------------------------------
    for(auto n : nums)
        REQUIRE( indexfn(nums, lambda(x, x == n)) == n - 1 );

    REQUIRE( indexfn(nums, lambda(x, x == 1000)) >= nums.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: filter
    //-------------------------------------------------------------------------
    res = filter(nums, lambda(x, x % 2 == 0));
    REQUIRE( res == Vec<int>{2, 4, 6, 8} );

    //-------------------------------------------------------------------------
    // TESTING METHOD: remove
    //-------------------------------------------------------------------------
    res = remove(nums, lambda(x, x % 2 == 1));
    REQUIRE( res == Vec<int>{2, 4, 6, 8} );

    //-------------------------------------------------------------------------
    // TESTING METHOD: unique
    //-------------------------------------------------------------------------
    res = unique(Vec<int>{ 3, 2, 3, 6, 8, 2, 5 });
    REQUIRE( res == Vec<int>{ 2, 3, 5, 6, 8 } );

    //-------------------------------------------------------------------------
    // TESTING METHOD: vectorize
    //-------------------------------------------------------------------------
    strs = vectorize(Map<String, int>{{"C", 1}, {"A", 2}, {"B", 3}}, lambda(x, x.first));

    REQUIRE( strs.size() == 3 );
    REQUIRE( contains(strs, "A") );
    REQUIRE( contains(strs, "B") );
    REQUIRE( contains(strs, "C") );

    res = vectorize(Map<String, int>{{"C", 1}, {"A", 2}, {"B", 3}}, lambda(x, x.second));

    REQUIRE( res.size() == 3 );
    REQUIRE( contains(res, 1) );
    REQUIRE( contains(res, 2) );
    REQUIRE( contains(res, 3) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: concatenate
    //-------------------------------------------------------------------------
    res = concatenate(Vec<int>{ 3, 2, 1 }, Vec<int>{ 5, 2, 3, 4 });
    REQUIRE( res == Vec<int>{ 3, 2, 1, 5, 2, 3, 4 } );

    res = concatenate(Vec<int>{ 5, 2, 3, 4 }, Vec<int>{ 3, 2, 1 } );
    REQUIRE( res == Vec<int>{ 5, 2, 3, 4, 3, 2, 1 } );

    res = concatenate(Vec<int>{ 1, 2, 3, 4 }, Vec<int>{ 7, 6, 5 } );
    REQUIRE( res == Vec<int>{ 1, 2, 3, 4, 7, 6, 5 } );

    res = concatenate(Vec<int>{ 1, 2, 3, 4 }, Vec<int>{} );
    REQUIRE( res == Vec<int>{ 1, 2, 3, 4 } );

    res = concatenate(Vec<int>{}, Vec<int>{1, 2, 3, 4} );
    REQUIRE( res == Vec<int>{ 1, 2, 3, 4 } );

    res = concatenate(Vec<int>{}, Vec<int>{} );
    REQUIRE( res == Vec<int>{} );

    //-------------------------------------------------------------------------
    // TESTING METHOD: merge
    //-------------------------------------------------------------------------
    res = merge(Vec<int>{ 3, 2, 1 }, Vec<int>{ 5, 2, 3, 4 });
    REQUIRE( res == Vec<int>{1, 2, 3, 4, 5} );

    res = merge(Vec<int>{ 5, 2, 3, 4 }, Vec<int>{ 3, 2, 1 } );
    REQUIRE( res == Vec<int>{1, 2, 3, 4, 5} );

    res = merge(Vec<int>{ 1, 2, 3, 4 }, Vec<int>{ 7, 6, 5 } );
    REQUIRE( res == Vec<int>{1, 2, 3, 4, 5, 6, 7} );

    res = merge(Vec<int>{ 1, 2, 3, 4 }, Vec<int>{} );
    REQUIRE( res == Vec<int>{1, 2, 3, 4} );

    res = merge(Vec<int>{}, Vec<int>{1, 2, 3, 4} );
    REQUIRE( res == Vec<int>{1, 2, 3, 4} );

    res = merge(Vec<int>{}, Vec<int>{} );
    REQUIRE( res == Vec<int>{} );

    //-------------------------------------------------------------------------
    // TESTING METHOD: contains
    //-------------------------------------------------------------------------
    REQUIRE(  contains(Vec<int>{1, 2, 3}, 3) );
    REQUIRE( !contains(Vec<int>{1, 2, 3}, 5) );

    REQUIRE(  contains(Strings{"1", "2", "3"}, "3") );
    REQUIRE( !contains(Strings{"1", "2", "3"}, "5") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: containsfn
    //-------------------------------------------------------------------------
    REQUIRE(  containsfn(Vec<int>{1, 2, 5}, lambda(x, x == 5)) );
    REQUIRE( !containsfn(Vec<int>{1, 2, 5}, lambda(x, x == 3)) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: contained
    //-------------------------------------------------------------------------
    REQUIRE(  contained(Vec<int>{},        Vec<int>{1, 2, 3, 4, 5}) );
    REQUIRE(  contained(Vec<int>{1},       Vec<int>{1, 2, 3, 4, 5}) );
    REQUIRE(  contained(Vec<int>{5, 3},    Vec<int>{1, 2, 3, 4, 5}) );
    REQUIRE(  contained(Vec<int>{2, 4, 5}, Vec<int>{1, 2, 3, 4, 5}) );
    REQUIRE( !contained(Vec<int>{6},       Vec<int>{1, 2, 3, 4, 5}) );
    REQUIRE( !contained(Vec<int>{2, 9, 5}, Vec<int>{1, 2, 3, 4, 5}) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: disjoint
    //-------------------------------------------------------------------------
    REQUIRE(  disjoint(Vec<int>{},        Vec<int>{1, 2, 3, 4, 5}) );
    REQUIRE(  disjoint(Vec<int>{6},       Vec<int>{1, 2, 3, 4, 5}) );
    REQUIRE(  disjoint(Vec<int>{6, 8},    Vec<int>{1, 2, 3, 4, 5}) );
    REQUIRE(  disjoint(Vec<int>{6, 8, 9}, Vec<int>{1, 2, 3, 4, 5}) );
    REQUIRE( !disjoint(Vec<int>{4},       Vec<int>{1, 2, 3, 4, 5}) );
    REQUIRE( !disjoint(Vec<int>{2, 1, 5}, Vec<int>{1, 2, 3, 4, 5}) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: identical
    //-------------------------------------------------------------------------
    REQUIRE(  identical(Vec<int>{},        Vec<int>{}) );
    REQUIRE(  identical(Vec<int>{6},       Vec<int>{6}) );
    REQUIRE(  identical(Vec<int>{6, 8},    Vec<int>{8, 6}) );
    REQUIRE(  identical(Vec<int>{6, 8, 9}, Vec<int>{9, 6, 8}) );
    REQUIRE( !identical(Vec<int>{4},       Vec<int>{5}) );
    REQUIRE( !identical(Vec<int>{2, 1, 5}, Vec<int>{2, 1, 6}) );
    REQUIRE( !identical(Vec<int>{2, 1, 5}, Vec<int>{2, 1, 4, 5}) );
}
