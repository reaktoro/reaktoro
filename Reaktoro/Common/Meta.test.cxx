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
#include <Reaktoro/Common/Meta.hpp>
#include <Reaktoro/Common/Types.hpp>
using namespace Reaktoro;

TEST_CASE("Testing Meta", "[Meta]")
{
    //-------------------------------------------------------------------------
    // TESTING METHOD: For
    //-------------------------------------------------------------------------
    Vec<int> vec(4);
    int sum = 0;

    For<0, 4>([&](auto i) constexpr { sum += i; });
    CHECK( sum == 6 );

    For<0, 4>([&](auto i) constexpr { vec[i] = i*i; });
    CHECK( vec == Vec<int>{0, 1, 4, 9} );

    For<4>([&](auto i) constexpr { vec[i] = i + 10; });
    CHECK( vec == Vec<int>{10, 11, 12, 13} );

    For<2, 6>([&](auto i) constexpr { vec[i-2] = i*i; });
    CHECK( vec == Vec<int>{4, 9, 16, 25} );
}
