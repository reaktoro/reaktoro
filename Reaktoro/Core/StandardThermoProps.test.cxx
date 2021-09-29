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
#include <Reaktoro/Core/StandardThermoProps.hpp>
using namespace Reaktoro;

TEST_CASE("Testing StandardThermoProps class", "[StandardThermoProps]")
{
    StandardThermoProps props;
    props.G0  = 10.0;
    props.H0  = 20.0;
    props.V0  = 30.0;
    props.VT0 = 40.0;
    props.VP0 = 50.0;
    props.Cp0 = 60.0;

    const auto T = 100.0;
    const auto P = 3000.0;

    const auto [G0, H0, V0, Cp0, VT0, VP0] = props;

    // NOTHING HERE YET TO TEST...
}
