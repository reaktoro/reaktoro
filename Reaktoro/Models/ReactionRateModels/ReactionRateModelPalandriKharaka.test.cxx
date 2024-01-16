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
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Reaction.hpp>
#include <Reaktoro/Extensions/Supcrt/SupcrtDatabase.hpp>
#include <Reaktoro/Models/ReactionRateModels/ReactionRateModelPalandriKharaka.hpp>
#include <Reaktoro/Utils/AqueousProps.hpp>
using namespace Reaktoro;

struct TestingConditions
{
    real temperature;
    real pressure;
    real dissolved_mineral;
    real pH;
    real pCO2;
    real surface_area;
};

TestingConditions conditions1 = { 25.0,  1.0, 1e-6, 7.0, 0.1, 0.4 };
TestingConditions conditions2 = { 50.0, 10.0, 1e-3, 4.0, 9.0, 1.2 };
TestingConditions conditions3 = { 75.0, 20.0, 1e-1, 9.0, 1.0, 3.2 };

TEST_CASE("Testing ReactionRateModelPalandriKharaka class", "[ReactionRateModelPalandriKharaka]")
{
    TestingConditions conditions = GENERATE(as<TestingConditions>{},
        conditions1,
        conditions2,
        conditions3
    );

    ReactionRateModelParamsPalandriKharaka params;
    params.mineral = "Calcite";
    params.mechanisms = {{
        { "Acid", -0.30, 14.4, 1.0, 1.0, {{"H+", "a", 1.0}} },
        { "Neutral", -5.81, 23.5, 1.0, 1.0 },
        { "Carbonate", -3.48, 35.4, 1.0, 1.0, {{"CO2", "P", 1.0}} }
    }};

    SupcrtDatabase db("supcrtbl");

    ChemicalSystem system(db,
        AqueousPhase("H2O(aq) H+ OH- Ca+2 HCO3- CO3-2 CO2(aq)"),
        GaseousPhase("CO2(g) N2(g)"),
        MineralPhase("Calcite"),
        GeneralReaction("Calcite").setRateModel(ReactionRateModelPalandriKharaka(params)),
        Surface("Calcite").withAreaModel([=](ChemicalProps const&) { return conditions.surface_area; })
    );

    ChemicalState state(system);
    state.temperature(conditions.temperature, "celsius");
    state.pressure(conditions.pressure, "bar");
    state.set("H2O(aq)", 1.0, "kg");
    state.set("H+", pow(10, -conditions.pH), "mol");
    state.set("OH-", pow(10, conditions.pH - 14), "mol");
    state.set("Ca+2", conditions.dissolved_mineral, "mol");
    state.set("CO3-2", conditions.dissolved_mineral, "mol");
    state.set("CO2(g)", conditions.pCO2/conditions.pressure, "mol");
    state.set("N2(g)", 1 - conditions.pCO2/conditions.pressure, "mol");
    state.set("Calcite", 1.0, "mol");

    ChemicalProps props(state);

    AqueousProps aprops(props);

    const auto T = state.temperature();
    const auto P = state.pressure();
    const auto R = universalGasConstant;

    const auto SA = conditions.surface_area;

    const auto T0 = 298.15;

    const auto k0_acid = pow(10.0, params.mechanisms[0].lgk);
    const auto k0_neutral = pow(10.0, params.mechanisms[1].lgk);
    const auto k0_carbonate = pow(10.0, params.mechanisms[2].lgk);

    const auto E_acid = params.mechanisms[0].E * 1e3; // from kJ to J
    const auto E_neutral = params.mechanisms[1].E * 1e3; // from kJ to J
    const auto E_carbonate = params.mechanisms[2].E * 1e3; // from kJ to J

    const auto k_acid = k0_acid * exp(-E_acid/R * (1/T - 1/T0));
    const auto k_neutral = k0_neutral * exp(-E_neutral/R * (1/T - 1/T0));
    const auto k_carbonate = k0_carbonate * exp(-E_carbonate/R * (1/T - 1/T0));

    const auto Omega = aprops.saturationRatio(params.mineral);

    const auto pOmega_acid = pow(Omega, params.mechanisms[0].p);
    const auto pOmega_neutral = pow(Omega, params.mechanisms[1].p);
    const auto pOmega_carbonate = pow(Omega, params.mechanisms[2].p);

    const auto qOmega_acid = pow(1 - pOmega_acid, params.mechanisms[0].q);
    const auto qOmega_neutral = pow(1 - pOmega_neutral, params.mechanisms[1].q);
    const auto qOmega_carbonate = pow(1 - pOmega_carbonate, params.mechanisms[2].q);

    const auto g_acid = pow(props.speciesActivity("H+"), params.mechanisms[0].catalysts[0].power);
    const auto g_carbonate = pow(props.speciesMoleFraction("CO2(g)") * P * 1e-5, params.mechanisms[2].catalysts[0].power);

    const auto rate_expected = SA * (
        k_acid * qOmega_acid * g_acid +
        k_neutral * qOmega_neutral +
        k_carbonate * qOmega_carbonate * g_carbonate
    );

    const auto rate_actual = system.reaction(0).rate(props);

    CHECK( rate_actual == Approx(rate_expected) );
}
