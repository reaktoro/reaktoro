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
#include <Reaktoro/Utils/MineralSurface.hpp>
using namespace Reaktoro;

TEST_CASE("Testing MineralSurface", "[MineralSurface]")
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
    state.scalePhaseAmount("Calcite", 2.0, "mol");
    state.scalePhaseMass("Magnesite", 2.0, "kg");
    state.scalePhaseVolume("Quartz" , 2.0, "m3");

    ChemicalProps props(state);

    WHEN("MineralSurface is set with a constant surface area model")
    {
        auto surface1 = MineralSurface("Calcite", 1.23, "m2");
        auto surface2 = MineralSurface("Magnesite", 1.23, "cm2");
        auto surface3 = MineralSurface("Quartz", 1.23, "mm2");

        CHECK( surface1.name() == "Calcite" );
        CHECK( surface2.name() == "Magnesite" );
        CHECK( surface3.name() == "Quartz" );

        CHECK( surface1.areaModel()(props).val() == Approx(1.23) );
        CHECK( surface2.areaModel()(props).val() == Approx(1.23e-4) );
        CHECK( surface3.areaModel()(props).val() == Approx(1.23e-6) );
    }

    WHEN("MineralSurface is set with a linear surface area model")
    {
        auto surface1 = MineralSurface("Calcite", 0.5, "cm2/mmol");
        auto surface2 = MineralSurface("Magnesite", 0.5, "cm2/g");
        auto surface3 = MineralSurface("Quartz", 0.5, "mm2/mm3");

        CHECK( surface1.name() == "Calcite" );
        CHECK( surface2.name() == "Magnesite" );
        CHECK( surface3.name() == "Quartz" );

        CHECK( surface1.areaModel()(props).val() == Approx(1.0e-1) );
        CHECK( surface2.areaModel()(props).val() == Approx(1.0e-1) );
        CHECK( surface3.areaModel()(props).val() == Approx(1.0e+3) );
    }

    WHEN("MineralSurface is set with a power law surface area model")
    {
        auto surface1 = MineralSurface("Calcite", 1000.0, "m2", 20.0e6, "umol", 3.0);
        auto surface2 = MineralSurface("Magnesite", 1000.0e4, "cm2", 20.0e3, "g", 3.0);
        auto surface3 = MineralSurface("Quartz", 1000.0e6, "mm2", 20.0e9, "mm3", 3.0);

        CHECK( surface1.name() == "Calcite" );
        CHECK( surface2.name() == "Magnesite" );
        CHECK( surface3.name() == "Quartz" );

        CHECK( surface1.areaModel()(props).val() == Approx(1.0) ); // 1000 * (2 / 20)**3
        CHECK( surface2.areaModel()(props).val() == Approx(1.0) ); // 1000 * (2 / 20)**3
        CHECK( surface3.areaModel()(props).val() == Approx(1.0) ); // 1000 * (2 / 20)**3
    }
}
