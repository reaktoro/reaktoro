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
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Math/Roots.hpp>
using namespace Reaktoro;

TEST_CASE("Testing cardano function in Roots module", "[Roots]")
{
    double b, c, d;

    auto check = [](double b, double c, double d)
    {
        const auto [x1, x2, x3] = cardano(b, c, d);
        auto checkzero = [&](auto x)
        {
            auto res = x*x*x + b*x*x + c*x + d;
            INFO(str("b = ", b, ", c = ", c, ", d = ", d));
            REQUIRE( res.real() == Approx(0.0).scale(1.0) );
            REQUIRE( res.imag() == Approx(0.0).scale(1.0) );
        };
        checkzero(x1);
        checkzero(x2);
        checkzero(x3);
    };

    Vec<double> coeffs = {-5.0, -3.0, -1.0, 0.0, 1.0, 3.0, 5.0};

    for(auto b : coeffs)
        for(auto c : coeffs)
            for(auto d : coeffs)
                check(b, c, d);
}

TEST_CASE("Testing newton function in Roots module", "[Roots]")
{
    std::function<std::tuple<double, double>(const double&)> f;

    //=========================================================================
    // Solve cos(x) = 0.0
    //=========================================================================
    {
        double x0 = 1.0;
        f = [](const double& x) { return std::make_tuple(std::cos(x), -std::sin(x)); };

        const auto [x, iters, success] = newton(f, x0, 1.0e-6, 100);

        const auto [fx, _] = f(x);

        INFO(str("newton failed with iters = ", iters, " and residual = ", fx));
        REQUIRE( success );
        REQUIRE( iters < 100 );
        REQUIRE( fx == Approx(0.0).scale(1.0) );
    }

    //=========================================================================
    // Solve exp(x) = 10.0
    //=========================================================================
    {
        double x0 = 1.0;
        f = [](const double& x) { return std::make_tuple(std::exp(x) - 10.0, std::exp(x)); };

        const auto [x, iters, success] = newton(f, x0, 1.0e-6, 100);

        const auto [fx, _] = f(x);

        INFO(str("newton failed with iters = ", iters, " and residual = ", fx));
        REQUIRE( success );
        REQUIRE( iters < 100 );
        REQUIRE( fx == Approx(0.0).scale(1.0) );
    }
}
