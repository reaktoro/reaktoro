// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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
#include <Reaktoro/Math/Roots.hpp>
using namespace Reaktoro;

TEST_CASE("Testing cardano function in Roots module", "[Roots]")
{
    double b, c, d;

    //=========================================================================
    // Solve x**3 - 5 = 0
    //=========================================================================
    {
        b = 0.0; c = 0.0; d = -5.0;
        const auto [x1, x2, x3] = cardano(b, c, d);
        REQUIRE( (x1*x1*x1).real() == Approx(5.0) );
        REQUIRE( (x1*x1*x1).imag() == Approx(0.0) );
    }

    //=========================================================================
    // Solve x**3 + 5 = 0
    //=========================================================================
    {
        b = 0.0; c = 0.0; d = 5.0;
        const auto [x1, x2, x3] = cardano(b, c, d);
        REQUIRE( (x1*x1*x1).real() == Approx(-5.0) );
        REQUIRE( (x1*x1*x1).imag() == Approx( 0.0) );
    }

    //=========================================================================
    // Solve x**3 + x**2 + x + 1 = 0
    //=========================================================================
    {
        b = 1.0; c = 1.0; d = 1.0;
        const auto [x1, x2, x3] = cardano(b, c, d);
        REQUIRE( x1.real() == Approx(-1.0) );
        REQUIRE( x1.imag() == Approx( 0.0) );
        REQUIRE( x2.real() == Approx( 0.0) );
        REQUIRE( x2.imag() == Approx(-1.0) );
        REQUIRE( x3.real() == Approx( 0.0) );
        REQUIRE( x3.imag() == Approx( 1.0) );
    }

    //=========================================================================
    // Solve x**3 - x**2 - x + 1 = 0
    //=========================================================================
    {
        b = -1.0; c = -1.0; d = 1.0;
        const auto [x1, x2, x3] = cardano(b, c, d);
        REQUIRE( x1.real() == Approx( 1.0) );
        REQUIRE( x1.imag() == Approx( 0.0) );
        REQUIRE( x2.real() == Approx(-1.0) );
        REQUIRE( x2.imag() == Approx( 0.0) );
        REQUIRE( x3.real() == Approx( 1.0) );
        REQUIRE( x3.imag() == Approx( 0.0) );
    }
}

TEST_CASE("Testing newton function in Roots module", "[Roots]")
{
    std::function<std::tuple<double, double>(const double&)> f;

    //=========================================================================
    // Solve cos(x) = 0.0
    //=========================================================================
    {
        double x0 = 0.0;
        f = [](const double& x) { return std::make_tuple(std::cos(x), -std::sin(x)); };

        const auto [x, iters, success] = newton(f, x0, 1.0e-6, 100);

        const auto [fx, _] = f(x);

        REQUIRE( success );
        REQUIRE( iters < 100 );
        REQUIRE( fx == Approx(0.0) );
    }

    //=========================================================================
    // Solve exp(x) = 10.0
    //=========================================================================
    {
        double x0 = 0.0;
        f = [](const double& x) { return std::make_tuple(std::exp(x) - 10.0, std::exp(x)); };

        const auto [x, iters, success] = newton(f, x0, 1.0e-6, 100);

        const auto [fx, _] = f(x);

        REQUIRE( success );
        REQUIRE( iters < 100 );
        REQUIRE( fx == Approx(0.0) );
    }
}
