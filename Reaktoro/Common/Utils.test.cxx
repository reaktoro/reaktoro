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
#include <Reaktoro/Common/Utils.hpp>
using namespace Reaktoro;

TEST_CASE("Testing Utils module", "[Utils]")
{
    // The possibilities for an array of species amounts
    auto n           = ArrayXr{{1.0, 2.0, 3.0, 4.0}};
    auto nzero       = ArrayXr{{0.0, 0.0, 0.0, 0.0}};
    auto nsingle     = ArrayXr{{2.0}};
    auto nsinglezero = ArrayXr{{0.0}};

    ArrayXr x;
    MatrixXd dxdn;

    //-------------------------------------------------------------------------
    // TESTING METHOD: moleFractions
    //-------------------------------------------------------------------------

    x = moleFractions(n);

    REQUIRE( x[0] == Approx(0.1) );
    REQUIRE( x[1] == Approx(0.2) );
    REQUIRE( x[2] == Approx(0.3) );
    REQUIRE( x[3] == Approx(0.4) );

    x = moleFractions(nzero);

    REQUIRE( x[0] == Approx(0.0) );
    REQUIRE( x[1] == Approx(0.0) );
    REQUIRE( x[2] == Approx(0.0) );
    REQUIRE( x[3] == Approx(0.0) );

    x = moleFractions(nsingle);

    REQUIRE( x[0] == Approx(1.0) );

    x = moleFractions(nsinglezero);

    REQUIRE( x[0] == Approx(1.0) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: moleFractionsJacobian
    //-------------------------------------------------------------------------
    dxdn = moleFractionsJacobian(n);

    INFO("dxdn = \n" << dxdn);
    REQUIRE( dxdn.isApprox(MatrixXd{{
        { 0.09, -0.01, -0.01, -0.01},
        {-0.02,  0.08, -0.02, -0.02},
        {-0.03, -0.03,  0.07, -0.03},
        {-0.04, -0.04, -0.04,  0.06},
    }}));

    dxdn = moleFractionsJacobian(nzero);

    INFO("dxdn = \n" << dxdn);
    REQUIRE( dxdn.isApprox(MatrixXd::Zero(4, 4)) );

    dxdn = moleFractionsJacobian(nsingle);

    INFO("dxdn = \n" << dxdn);
    REQUIRE( dxdn.isApprox(MatrixXd::Zero(1, 1)) );

    dxdn = moleFractionsJacobian(nsinglezero);

    INFO("dxdn = \n" << dxdn);
    REQUIRE( dxdn.isApprox(MatrixXd::Zero(1, 1)) );
}
