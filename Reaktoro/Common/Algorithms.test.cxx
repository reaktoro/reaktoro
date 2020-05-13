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
        REQUIRE( indexfn(nums, RKT_LAMBDA(x, x == n)) == n - 1 );

    REQUIRE( indexfn(nums, RKT_LAMBDA(x, x == 1000)) >= nums.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: filter
    //-------------------------------------------------------------------------
    res = filter(nums, RKT_LAMBDA(x, x % 2 == 0));
    REQUIRE( res == Vec<int>{2, 4, 6, 8} );

    //-------------------------------------------------------------------------
    // TESTING METHOD: remove
    //-------------------------------------------------------------------------
    res = remove(nums, RKT_LAMBDA(x, x % 2 == 1));
    REQUIRE( res == Vec<int>{2, 4, 6, 8} );

    //-------------------------------------------------------------------------
    // TESTING METHOD: unique
    //-------------------------------------------------------------------------
    res = unique(Vec<int>{ 3, 2, 3, 6, 8, 2, 5 });
    REQUIRE( res == Vec<int>{ 2, 3, 5, 6, 8 } );

    //-------------------------------------------------------------------------
    // TESTING METHOD: extract
    //-------------------------------------------------------------------------
    REQUIRE( extract(Strings{"A", "B", "C", "D"}, Vec<int>{}) == Strings{} );
    REQUIRE( extract(Strings{"A", "B", "C", "D"}, Vec<int>{0}) == Strings{"A"} );
    REQUIRE( extract(Strings{"A", "B", "C", "D"}, Vec<int>{0, 2}) == Strings{"A", "C"} );
    REQUIRE( extract(Strings{"A", "B", "C", "D"}, Vec<int>{0, 2, 1}) == Strings{"A", "C", "B"} );
    REQUIRE( extract(Strings{"A", "B", "C", "D"}, Vec<int>{0, 2, 1, 2}) == Strings{"A", "C", "B", "C"} );

    //-------------------------------------------------------------------------
    // TESTING METHOD: vectorize
    //-------------------------------------------------------------------------
    strs = vectorize(Map<String, int>{{"C", 1}, {"A", 2}, {"B", 3}}, RKT_LAMBDA(x, x.first));

    REQUIRE( strs.size() == 3 );
    REQUIRE( contains(strs, "A") );
    REQUIRE( contains(strs, "B") );
    REQUIRE( contains(strs, "C") );

    res = vectorize(Map<String, int>{{"C", 1}, {"A", 2}, {"B", 3}}, RKT_LAMBDA(x, x.second));

    REQUIRE( res.size() == 3 );
    REQUIRE( contains(res, 1) );
    REQUIRE( contains(res, 2) );
    REQUIRE( contains(res, 3) );

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
    REQUIRE(  containsfn(Vec<int>{1, 2, 5}, RKT_LAMBDA(x, x == 5)) );
    REQUIRE( !containsfn(Vec<int>{1, 2, 5}, RKT_LAMBDA(x, x == 3)) );

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
    // TESTING METHOD: intersect
    //-------------------------------------------------------------------------
    REQUIRE( intersect(Vec<int>{}, Vec<int>{}) == Vec<int>{} );
    REQUIRE( intersect(Vec<int>{1, 2, 3}, Vec<int>{}) == Vec<int>{} );
    REQUIRE( intersect(Vec<int>{}, Vec<int>{1, 2, 3}) == Vec<int>{} );
    REQUIRE( intersect(Vec<int>{1, 2, 3}, Vec<int>{3, 2, 1}) == Vec<int>{1, 2, 3} );
    REQUIRE( intersect(Vec<int>{3, 2, 1}, Vec<int>{1, 2, 3}) == Vec<int>{3, 2, 1} );
    REQUIRE( intersect(Vec<int>{2, 1}, Vec<int>{1, 2, 3}) == Vec<int>{2, 1} );
    REQUIRE( intersect(Vec<int>{2, 1}, Vec<int>{4, 5, 6}) == Vec<int>{} );

    //-------------------------------------------------------------------------
    // TESTING METHOD: difference
    //-------------------------------------------------------------------------
    REQUIRE( difference(Vec<int>{}, Vec<int>{}) == Vec<int>{} );
    REQUIRE( difference(Vec<int>{1, 2, 3}, Vec<int>{}) == Vec<int>{1, 2, 3} );
    REQUIRE( difference(Vec<int>{}, Vec<int>{1, 2, 3}) == Vec<int>{} );
    REQUIRE( difference(Vec<int>{1, 2, 3}, Vec<int>{3, 2, 1}) == Vec<int>{} );
    REQUIRE( difference(Vec<int>{3, 2, 1}, Vec<int>{1, 2, 3}) == Vec<int>{} );
    REQUIRE( difference(Vec<int>{2, 1}, Vec<int>{1, 2, 3}) == Vec<int>{} );
    REQUIRE( difference(Vec<int>{2, 1}, Vec<int>{4, 5, 6}) == Vec<int>{2, 1} );

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

    //-------------------------------------------------------------------------
    // TESTING METHOD: range
    //-------------------------------------------------------------------------
    REQUIRE( range(1, 1, 2) == Vec<int>{} );
    REQUIRE( range(1, 5, 2) == Vec<int>{1, 3} );
    REQUIRE( range(1, 6, 2) == Vec<int>{1, 3, 5} );
    REQUIRE( range(1, 4) == Vec<int>{1, 2, 3} );
    REQUIRE( range(0) == Vec<int>{} );
    REQUIRE( range(3) == Vec<int>{0, 1, 2} );

    //-------------------------------------------------------------------------
    // TESTING METHOD: oneof
    //-------------------------------------------------------------------------
    REQUIRE( oneof(1, 1) );
    REQUIRE( oneof(1, 3, 4, 5, 1, 7) );

    REQUIRE_FALSE( oneof(1, 2) );
    REQUIRE_FALSE( oneof(1, 4, 5) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: sum
    //-------------------------------------------------------------------------
    REQUIRE( sum(2, 5, [](int i) { return i*i; }) == Approx(2*2 + 3*3 + 4*4) );
    REQUIRE( sum(3, [](int i) { return i*i; }) == Approx(0*0 + 1*1 + 2*2) );
}
