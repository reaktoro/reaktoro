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
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Extensions/Supcrt/SupcrtDatabase.hpp>
#include <Reaktoro/Models/SurfaceAreaModels/SurfaceAreaModelLinear.hpp>
using namespace Reaktoro;

TEST_CASE("Testing SurfaceAreaModelLinear", "[SurfaceAreaModelLinear]")
{
    SupcrtDatabase db("supcrtbl");

    ChemicalSystem system(db,
        AqueousPhase("H2O(aq) H+ OH- Ca+2 HCO3- CO2(aq) CO3-2"),
        MineralPhase("Calcite"),
        MineralPhase("Magnesite"),
        MineralPhase("Quartz")
    );

    ChemicalState state(system);

    state.set("H2O(aq)", 1.0, "kg");
    state.scalePhaseAmount("Calcite", 10.0, "mol");
    state.scalePhaseMass("Magnesite", 10.0, "kg");
    state.scalePhaseVolume("Quartz" , 10.0, "m3");

    ChemicalProps props(state);

    auto model1 = SurfaceAreaModelLinearMolar("Calcite", 0.1);
    auto model2 = SurfaceAreaModelLinearSpecific("Magnesite", 0.1);
    auto model3 = SurfaceAreaModelLinearSpecific("Quartz", 0.1);

    CHECK( model1(props) == Approx(1.0) );
    CHECK( model2(props) == Approx(1.0) );
    CHECK( model3(props) == Approx(1.0) );
}
