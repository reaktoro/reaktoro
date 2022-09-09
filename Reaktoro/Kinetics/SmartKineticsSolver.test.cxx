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
#include <Reaktoro/Equilibrium/EquilibriumConditions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumRestrictions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>
#include <Reaktoro/Equilibrium/EquilibriumUtils.hpp>
#include <Reaktoro/Extensions/Supcrt/SupcrtDatabase.hpp>
#include <Reaktoro/Kinetics/SmartKineticsOptions.hpp>
#include <Reaktoro/Kinetics/SmartKineticsResult.hpp>
#include <Reaktoro/Kinetics/SmartKineticsSolver.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelDavies.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelPitzerHMW.hpp>

using namespace Reaktoro;

/// Return the largest relative difference between two arrays `actual` and `expected`.
auto largestRelativeDifference(ArrayXr const& actual, ArrayXr const& expected) -> double
{
    return ((actual - expected)/expected).abs().maxCoeff();
}

TEST_CASE("Testing SmartKineticsSolver", "[SmartKineticsSolver]")
{
    WHEN("temperature and pressure are given - calcite and water")
    {
        SupcrtDatabase db("supcrtbl");

        AqueousPhase solution("H2O(aq) H+ OH- Ca+2 HCO3- CO3-2 CO2(aq)");
        solution.setActivityModel(ActivityModelPitzerHMW());

        MineralPhase calcite("Calcite");

        ChemicalSystem system(db, solution, calcite);

        SmartKineticsResult result;
        SmartKineticsSolver solver(system);

        ChemicalState exactstate(system);
        KineticsSolver exactsolver(system);

        //-------------------------------------------------------------------------------------------------------------
        // CREATE AN INITIAL CHEMICAL STATE AND EQUILIBRATE IT - THIS IS THE FIRST LEARNING OPERATION
        //-------------------------------------------------------------------------------------------------------------

        ChemicalState state(system);
        state.temperature(25.0, "celsius");
        state.pressure(1.0, "bar");
        state.set("H2O(aq)", 1.0, "kg");
        state.set("Calcite", 1.0, "mol");

        result = solver.solve(state);

        CHECK( result.succeeded() );
        CHECK( result.learned() );
        CHECK( result.iterations() == 17 );

        //-------------------------------------------------------------------------------------------------------------
        // CHANGE THE CHEMICAL STATE SLIGHTLY AND CHECK SMART PREDICTION SUCCEEDED
        //-------------------------------------------------------------------------------------------------------------

        state.temperature(28.0, "celsius");
        state.pressure(2.0, "bar");
        state.set("H2O(aq)", 1.1, "kg");
        state.set("Calcite", 1.1, "mol");

        exactstate = state;
        equilibrate(exactstate); // compute the equilibrium state exactly with conventional algorithm

        result = solver.solve(state);

        CHECK( result.succeeded() );
        CHECK( result.predicted() );
        CHECK( result.iterations() == 0 );

        CHECK( largestRelativeDifference(state.speciesAmounts(), exactstate.speciesAmounts()) == Approx(0.0285944) ); // ~2.8% max relative difference

        //-------------------------------------------------------------------------------------------------------------
        // CHANGE THE CHEMICAL STATE MORE STRONGLY AND CHECK A LEARNING OPERATION WAS NEEEDED
        //-------------------------------------------------------------------------------------------------------------

        state.temperature(50.0, "celsius");
        state.pressure(10.0, "bar");
        state.set("H2O(aq)", 2.0, "kg");
        state.set("Calcite", 2.0, "mol");

        exactstate = state;
        exactsolver.solve(exactstate);

        result = solver.solve(state);

        CHECK( result.succeeded() );
        CHECK( result.learned() );
        CHECK( result.iterations() == 6 );
    }
}
