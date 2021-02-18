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
#include <Reaktoro/Core/Param.hpp>
using namespace Reaktoro;

TEST_CASE("Testing Param class", "[Param]")
{
    Param x;

    CHECK( x.value() == 0.0 );
    CHECK( x.name() == "" );
    CHECK( x.lowerbound() < 0.0 );
    CHECK( x.upperbound() > 0.0 );
    CHECK( std::isinf(x.lowerbound()) );
    CHECK( std::isinf(x.upperbound()) );
    CHECK( x.isconst() == false );

    x = Param(1.0);

    CHECK( x.value() == 1.0 );
    CHECK( x.name() == "" );
    CHECK( x.lowerbound() < 0.0 );
    CHECK( x.upperbound() > 0.0 );
    CHECK( std::isinf(x.lowerbound()) );
    CHECK( std::isinf(x.upperbound()) );
    CHECK( x.isconst() == false );

    x = Param().value(3.0).name("x").lowerbound(1.0).upperbound(7.0).isconst(true);

    CHECK( x.value() == 3.0 );
    CHECK( x.name() == "x" );
    CHECK( x.lowerbound() == 1.0 );
    CHECK( x.upperbound() == 7.0 );
    CHECK( x.isconst() == true );

    x.lowerbound(0.0).name("y").upperbound(17.0).value(9.0).isconst(false);

    CHECK( x.value() == 9.0 );
    CHECK( x.name() == "y" );
    CHECK( x.lowerbound() == 0.0 );
    CHECK( x.upperbound() == 17.0 );
    CHECK( x.isconst() == false );

    Param z = x; // z points to x

    CHECK( z.value() == 9.0 );
    CHECK( z.name() == "y" );
    CHECK( z.lowerbound() == 0.0 );
    CHECK( z.upperbound() == 17.0 );
    CHECK( z.isconst() == false );

    z.value(10.0).name("zzz"); // changing Param z here actually changes underlying Param x

    CHECK( x.value() == 10.0 ); // x.value() changed with z.value(10.0)
    CHECK( x.name() == "zzz" );   // x.name() changed with z.name("z")

    x = Param::Constant(11.0);

    CHECK( x.value() == 11.0 );
    CHECK( x.name() == "" );
    CHECK( x.lowerbound() < 0.0 );
    CHECK( x.upperbound() > 0.0 );
    CHECK( std::isinf(x.lowerbound()) );
    CHECK( std::isinf(x.upperbound()) );
    CHECK( x.isconst() == true );
}
