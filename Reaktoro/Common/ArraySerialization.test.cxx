// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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
#include <Reaktoro/Common/ArraySerialization.hpp>
using namespace Reaktoro;

TEST_CASE("Testing ArraySerialization", "[ArraySerialization]")
{
    ArrayXd array;

    ArraySerialization::resize(array, 1.0, 2.0);
    CHECK( array.size() == 2 );

    ArraySerialization::serialize(array, 1.0, 2.0);
    CHECK( array[0] == 1.0 );
    CHECK( array[1] == 2.0 );

    ArraySerialization::resize(array, 1.0, 2.0, ArrayXd{{3.0, 4.0, 5.0}}, 6.0, 7.0);
    CHECK( array.size() == 7 );

    ArraySerialization::serialize(array, 1.0, 2.0, ArrayXd{{3.0, 4.0, 5.0}}, 6.0, 7.0);
    CHECK( array[0] == 1.0 );
    CHECK( array[1] == 2.0 );
    CHECK( array[2] == 3.0 );
    CHECK( array[3] == 4.0 );
    CHECK( array[4] == 5.0 );
    CHECK( array[5] == 6.0 );
    CHECK( array[6] == 7.0 );

    double a;
    double b;
    VectorXd c(3);
    VectorXd d(2);

    ArraySerialization::deserialize(array, a, b, c, d);
    CHECK( a == 1.0 );
    CHECK( b == 2.0 );
    CHECK( c == VectorXd{{3.0, 4.0, 5.0}} );
    CHECK( d == VectorXd{{6.0, 7.0}} );
}
