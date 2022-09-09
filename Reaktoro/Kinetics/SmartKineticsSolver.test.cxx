// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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

// C++ includes
#include <iostream>

// Catch includes
#include <catch2/catch.hpp>

// Reaktoro includes
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Params.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConditions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumRestrictions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>
#include <Reaktoro/Extensions/Supcrt/SupcrtDatabase.hpp>
#include <Reaktoro/Kinetics/KineticsResult.hpp>
#include <Reaktoro/Kinetics/KineticsSolver.hpp>
#include <Reaktoro/Kinetics/SmartKineticsOptions.hpp>
#include <Reaktoro/Kinetics/SmartKineticsResult.hpp>
#include <Reaktoro/Kinetics/SmartKineticsSolver.hpp>
#include <Reaktoro/Math/MathUtils.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelDavies.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelPitzerHMW.hpp>
#include <Reaktoro/Models/ReactionRateModels/ReactionRateModelPalandriKharaka.hpp>
#include <Reaktoro/Models/ReactionRateModels/Support/MineralReactions.hpp>
using namespace Reaktoro;

TEST_CASE("Testing SmartKineticsSolver", "[SmartKineticsSolver]")
{
    WHEN("temperature and pressure are given - calcite and water")
    {
        SupcrtDatabase db("supcrtbl");

        AqueousPhase solution("H2O(aq) H+ OH- Ca+2 HCO3- CO3-2 CO2(aq)");
        solution.set(ActivityModelDavies());

        MineralPhase calcite("Calcite");

        Phases phases(db);
        phases.add(solution);
        phases.add(calcite);

        Params params = Params::embedded("PalandriKharaka.yaml");

        MineralReactions reactions("Calcite");
        reactions.setRateModel(ReactionRateModelPalandriKharaka(params));

        ChemicalSystem system(phases, reactions);

        ChemicalState state(system);
        ChemicalState exactstate(system);

        KineticsSolver exactsolver(system);
        SmartKineticsSolver solver(system);

        SmartKineticsResult result;

        real dt; // time step in seconds

        //-------------------------------------------------------------------------------------------------------------
        // CREATE AN INITIAL CHEMICAL STATE AND REACT IT - THIS IS THE FIRST LEARNING OPERATION
        //-------------------------------------------------------------------------------------------------------------

        dt = 0.1;

        state = ChemicalState(system);
        state.temperature(25.0, "celsius");
        state.pressure(1.0, "bar");
        state.set("H2O(aq)", 1.0, "kg");
        state.set("Calcite", 1.0, "mol");
        state.surfaceArea("Calcite", 1.0, "m2");

        result = solver.solve(state, dt);

        CHECK( result.succeeded() );
        CHECK( result.learned() );
        CHECK( result.iterations() == 15 );

        //-------------------------------------------------------------------------------------------------------------
        // CHANGE THE INITIAL CHEMICAL STATE SLIGHTLY AND CHECK SMART PREDICTION SUCCEEDED
        //-------------------------------------------------------------------------------------------------------------

        dt = 0.12;

        state = ChemicalState(system);
        state.temperature(30.0, "celsius");
        state.pressure(2.0, "bar");
        state.set("H2O(aq)", 1.1, "kg");
        state.set("Calcite", 1.1, "mol");
        state.surfaceArea("Calcite", 1.1, "m2");

        exactstate = state;
        exactsolver.solve(exactstate, dt); // compute the reacted state exactly with conventional algorithm

        result = solver.solve(state, dt);

        CHECK( result.succeeded() );
        CHECK( result.predicted() );
        CHECK( result.iterations() == 0 );

        CHECK( largestRelativeDifference(state.speciesAmounts(), exactstate.speciesAmounts()) == Approx(0.0577676443) ); // ~5.8% max relative difference
        CHECK( largestRelativeDifferenceLogScale(state.speciesAmounts(), exactstate.speciesAmounts()) == Approx(0.0051315941) ); // ~0.5% max relative difference

        //-------------------------------------------------------------------------------------------------------------
        // CHANGE THE INITIAL CHEMICAL STATE MORE STRONGLY AND CHECK A LEARNING OPERATION WAS NEEEDED
        //-------------------------------------------------------------------------------------------------------------

        dt = 0.15;

        state = ChemicalState(system);
        state.temperature(50.0, "celsius");
        state.pressure(10.0, "bar");
        state.set("H2O(aq)", 2.0, "kg");
        state.set("Calcite", 2.0, "mol");
        state.surfaceArea("Calcite", 1.0, "m2");

        exactstate = state;
        exactsolver.solve(exactstate, dt);

        result = solver.solve(state, dt);

        CHECK( result.succeeded() );
        CHECK( result.learned() );
        CHECK( result.iterations() == 16 );
    }
}
