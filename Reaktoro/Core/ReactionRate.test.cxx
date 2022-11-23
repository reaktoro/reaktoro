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
#include <Reaktoro/Core/ReactionRate.hpp>
using namespace Reaktoro;

TEST_CASE("Testing ReactionRate", "[ReactionRate]")
{
    auto rate = ReactionRate();

    rate = ReactionRate(1.23);
    CHECK( rate.value() == 1.23 );
    CHECK( rate.onEquationMode() == false );

    rate = ReactionRate(real(2.34));
    CHECK( rate.value() == 2.34 );
    CHECK( rate.onEquationMode() == false );

    rate = ReactionRate::enforce(123.4);
    CHECK( rate.value() == 123.4 );
    CHECK( rate.onEquationMode() == true );

    rate = 1.0;
    CHECK( rate.value() == 1.0 );
    CHECK( rate.onEquationMode() == true );

    rate = real(2.0);
    CHECK( rate.value() == 2.0 );
    CHECK( rate.onEquationMode() == true );

    SECTION("Testing arithmetic assign-operators with double and real scalars")
    {
        auto seed = GENERATE(1.0, real(1.0));

        rate = ReactionRate(2.0 * seed);

        rate += (1.0 * seed);
        CHECK( rate.value() == Approx(3.0) );

        rate -= (1.0 * seed);
        CHECK( rate.value() == Approx(2.0) );

        rate *= (3.0 * seed);
        CHECK( rate.value() == Approx(6.0) );

        rate /= (3.0 * seed);
        CHECK( rate.value() == Approx(2.0) );
    }

    SECTION("Testing positive and negative arithmetic operators")
    {
        rate = ReactionRate(2.0);

        rate = +rate;
        CHECK( rate.value() == Approx(2.0) );

        rate = -rate;
        CHECK( rate.value() == Approx(-2.0) );
    }

    SECTION("Testing arithmetic operators with double and real scalars on the right")
    {
        auto seed = GENERATE(1.0, real(1.0));

        rate = ReactionRate(2.0 * seed);

        rate = rate + (1.0 * seed);
        CHECK( rate.value() == Approx(3.0) );

        rate = rate - (1.0 * seed);
        CHECK( rate.value() == Approx(2.0) );

        rate = rate * (3.0 * seed);
        CHECK( rate.value() == Approx(6.0) );

        rate = rate / (3.0 * seed);
        CHECK( rate.value() == Approx(2.0) );
    }

    SECTION("Testing arithmetic operators with double and real scalars on the left")
    {
        auto seed = GENERATE(1.0, real(1.0));

        rate = ReactionRate(2.0 * seed);

        rate = (1.0 * seed) + rate;
        CHECK( rate.value() == Approx(3.0) );

        rate = (1.0 * seed) - rate;
        CHECK( rate.value() == Approx(2.0) );

        rate = (3.0 * seed) * rate;
        CHECK( rate.value() == Approx(6.0) );

        rate = (3.0 * seed) / rate;
        CHECK( rate.value() == Approx(2.0) );
    }
}
