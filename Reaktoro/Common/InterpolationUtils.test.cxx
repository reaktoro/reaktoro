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
#include <Reaktoro/Common/InterpolationUtils.hpp>
using namespace Reaktoro;

TEST_CASE("Testing InterpolationUtils module", "[InterpolationUtils]")
{
    SECTION("Checking method interpolateLinear")
    {
        const auto a = 2.3; // a in L(x) = a + bx
        const auto b = 6.7; // b in L(x) = a + bx

        const auto L = [&](auto x) { return a + b*x; };

        const auto x0 = 5.7;
        const auto x1 = 9.2;

        const auto y0 = L(x0);
        const auto y1 = L(x1);

        for(auto x : Vec<double>{-12, -8.0, -4.0, 0.0, 4.0, 8.0, 12})
            CHECK( interpolateLinear(x, x0, x1, y0, y1) == Approx(L(x)) );
    }

    SECTION("Checking method interpolateQuadratic")
    {
        const auto a = 2.3; // a in L(x) = a + bx + cx^2
        const auto b = 6.7; // b in L(x) = a + bx + cx^2
        const auto c = 8.4; // x in L(x) = a + bx + cx^2

        const auto Q = [&](auto x) { return a + b*x + c*x*x; };

        const auto x0 = 5.7;
        const auto x1 = 7.2;
        const auto x2 = 9.4;

        const auto y0 = Q(x0);
        const auto y1 = Q(x1);
        const auto y2 = Q(x2);

        for(auto x : Vec<double>{-12, -8.0, -4.0, 0.0, 4.0, 8.0, 12})
            CHECK( interpolateQuadratic(x, x0, x1, x2, y0, y1, y2) == Approx(Q(x)) );
    }
}
