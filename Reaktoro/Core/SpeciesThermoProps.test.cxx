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
#include <Reaktoro/Core/SpeciesThermoProps.hpp>
#include <Reaktoro/Core/StandardThermoProps.hpp>
using namespace Reaktoro;

TEST_CASE("Testing SpeciesThermoProps class", "[SpeciesThermoProps]")
{
    StandardThermoProps sprops;
    sprops.G0  = 10.0;
    sprops.H0  = 20.0;
    sprops.V0  = 30.0;
    sprops.VT0 = 40.0;
    sprops.VP0 = 50.0;
    sprops.Cp0 = 60.0;

    const auto T = 100.0;
    const auto P = 3000.0;

    const auto [G0, H0, V0, Cp0, VT0, VP0] = sprops;

    SpeciesThermoProps props(T, P, sprops);

    CHECK( props.T   == Approx(T)                   );
    CHECK( props.P   == Approx(P)                   );
    CHECK( props.G0  == Approx(G0)                  );
    CHECK( props.H0  == Approx(H0)                  );
    CHECK( props.V0  == Approx(V0)                  );
    CHECK( props.VT0 == Approx(VT0)                 );
    CHECK( props.VP0 == Approx(VP0)                 );
    CHECK( props.Cp0 == Approx(Cp0)                 );
    CHECK( props.Cv0 == Approx(Cp0 + T*VT0*VT0/VP0) );
    CHECK( props.U0  == Approx(H0 - P*V0)           );
    CHECK( props.S0  == Approx((H0 - G0)/T)         );
    CHECK( props.A0  == Approx(G0 - P*V0)           );
}
