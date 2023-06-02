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
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/Phases.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelPitzer.hpp>
#include <Reaktoro/Extensions/Phreeqc/PhreeqcDatabase.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelPhreeqcIonicStrengthPressureCorrection.hpp>
using namespace Reaktoro;

TEST_CASE("Testing ActivityModelPhreeqcIonicStrengthPressureCorrection", "[ActivityModelPhreeqcIonicStrengthPressureCorrection]")
{
    PhreeqcDatabase db("phreeqc.dat");

    AqueousPhase solution(speciate("Na Cl C"));
    solution.set(chain(
        ActivityModelPitzer(),
        ActivityModelPhreeqcIonicStrengthPressureCorrection()
    ));

    ChemicalSystem system(db, solution);

    ChemicalState state(system);
    state.temperature(50.0, "°C");
    state.pressure(100.0, "bar");
    state.set("H2O", 1.0, "kg");
    state.set("Na+", 4.0, "mol");
    state.set("Cl-", 4.0, "mol");
    state.set("Cl-", 4.0, "mol");

    state.props().update(state);

    auto ln_a = state.props().speciesActivitiesLn();

    CHECK( ln_a[0]  == Approx(-40.17420) );
    CHECK( ln_a[1]  == Approx(-36.14260) );
    CHECK( ln_a[2]  == Approx(-0.162252) );
    CHECK( ln_a[3]  == Approx(-36.20130) );
    CHECK( ln_a[4]  == Approx(-36.84130) );
    CHECK( ln_a[5]  == Approx(-37.89070) );
    CHECK( ln_a[6]  == Approx(-36.84130) );
    CHECK( ln_a[7]  == Approx(  1.17181) );
    CHECK( ln_a[8]  == Approx(-36.84130) );
    CHECK( ln_a[9]  == Approx(  1.17163) );
    CHECK( ln_a[10] == Approx(-37.94700) );
    CHECK( ln_a[11] == Approx(-36.84130) );
    CHECK( ln_a[12] == Approx(-37.44610) );
    CHECK( ln_a[13] == Approx(-36.84130) );
    CHECK( ln_a[14] == Approx(-36.84130) );
}
