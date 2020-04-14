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
using namespace Reaktoro;

TEST_CASE("Testing Algorithms", "[Algorithms]")
{
    std::vector<int> nums = { 1, 2, 3, 4, 5, 6, 7, 8, 9};
    std::vector<int> res;

    // Test the function index(container, value)
    for(auto n : nums)
        REQUIRE( index(nums, n) == n - 1);

    REQUIRE( index(nums, 1000) == -1 );

    // Test the function indexfn(container, pred)
    for(auto n : nums)
        REQUIRE( indexfn(nums, [=](auto x) { return x == n; } ) == n - 1);

    REQUIRE( indexfn(nums, [=](auto x) { return x == 1000; } ) == -1 );

    // Test the filter function
    res = filter(nums, [](auto x) { return x % 2 == 0; });

    REQUIRE( res.size() == 4 );
    REQUIRE( res[0] == 2 );
    REQUIRE( res[1] == 4 );
    REQUIRE( res[2] == 6 );
    REQUIRE( res[3] == 8 );

    // Test the remove function
    res = remove(nums, [](auto x) { return x % 2 == 1; });

    REQUIRE( res.size() == 4 );
    REQUIRE( res[0] == 2 );
    REQUIRE( res[1] == 4 );
    REQUIRE( res[2] == 6 );
    REQUIRE( res[3] == 8 );

    // Test the unique function
    std::vector<int> notunique = { 3, 2, 3, 6, 8, 2, 5 };

    res = unique(notunique);

    REQUIRE( res.size() == 5 );
    REQUIRE( res[0] == 2 );
    REQUIRE( res[1] == 3 );
    REQUIRE( res[2] == 5 );
    REQUIRE( res[3] == 6 );
    REQUIRE( res[4] == 8 );

    // Test the vectorize function
    std::unordered_map<std::string, int> map { {"C", 1}, {"A", 2}, {"B", 3} };

    res = vectorize(map, lambda(x, x.second));

    REQUIRE( res.size() == 3 );
    REQUIRE( contains(res, 1) );
    REQUIRE( contains(res, 2) );
    REQUIRE( contains(res, 3) );

    // Test the merge function
    std::vector<int> vec1 = { 3, 2, 1 };
    std::vector<int> vec2 = { 5, 2, 3, 4 };

    res = merge(vec1, vec2);

    REQUIRE( res.size() == 5 );
    REQUIRE( res[0] == 1 );
    REQUIRE( res[1] == 2 );
    REQUIRE( res[2] == 3 );
    REQUIRE( res[3] == 4 );
    REQUIRE( res[4] == 5 );

    // Test the contains function
    REQUIRE( contains(nums, 3) );
    REQUIRE_FALSE( contains(nums, 10) );

    // Test the containsfn function
    REQUIRE( containsfn(nums, [](auto x) { return x == 5; }) );
    REQUIRE_FALSE( containsfn(nums, [](auto x) { return x == 13; }) );

    // Test the contained function
    std::vector<int> a = {1, 2, 3, 4, 5};
    std::vector<int> b = {2, 4};
    std::vector<int> c = {4, 6};
    std::vector<int> d = {5, 7};
    std::vector<int> e = {};

    REQUIRE( contained(a, a) );
    REQUIRE( contained(b, a) );
    REQUIRE_FALSE( contained(c, a) );

    REQUIRE( disjoint(b, d) );
    REQUIRE( disjoint(b, e) );
    REQUIRE( disjoint(e, e) );
    REQUIRE_FALSE( disjoint(b, b) );
    REQUIRE_FALSE( disjoint(b, c) );
}
