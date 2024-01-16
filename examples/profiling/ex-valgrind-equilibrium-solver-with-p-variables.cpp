// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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

// -----------------------------------------------------------------------------
// 👏 Acknowledgements 👏
// -----------------------------------------------------------------------------
// This example was originally authored by:
//   • Allan Leal (30 September 2022)
// -----------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------
// Compile Reaktoro in Debug mode with optimization (e.g., -O2).
// Execute the command below:
//
// valgrind --tool=callgrind --collect-atstart=no examples/profiling/ex-valgrind-equilibrium-solver
//--------------------------------------------------------------------------------------------------

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

#include <valgrind/callgrind.h>

int main(int argc, char const *argv[])
{
    PhreeqcDatabase db("llnl.dat");

    AqueousPhase solution(speciate("C"));

    GaseousPhase gases("CO2(g) H2O(g)");

    ChemicalSystem system(db, solution, gases);

    ChemicalState waterState(system);
    waterState.set("H2O", 1.0, "kg");
    waterState.set("CO2", 50, "mg");

    EquilibriumSpecs specs(system);
    specs.temperature();
    specs.pressure();
    specs.pH();

    EquilibriumSolver solver(specs);

    EquilibriumConditions conditions(specs);
    conditions.temperature(7.2, "celsius");
    conditions.pressure(1.0, "atm");
    conditions.pH(8.0);

    solver.solve(waterState, conditions);

    EquilibriumSpecs specs2(system);
    specs2.temperature();
    specs2.phaseAmount("GaseousPhase");

    EquilibriumSolver solver2(specs2);

    EquilibriumConditions conditions2(specs2);
    conditions2.phaseAmount("GaseousPhase", 1, "umol");
    conditions2.setLowerBoundPressure(0.0005, "bar");
    conditions2.setUpperBoundPressure(200.0, "bar");

    conditions2.temperature(10.0, "celsius");

    ChemicalState state2(waterState);
    state2.add("CO2", 0.02, "kg");

    CALLGRIND_TOGGLE_COLLECT; // turn on
    auto res = solver2.solve(state2, conditions2);
    CALLGRIND_TOGGLE_COLLECT; // turn off

    errorifnot(res.succeeded(), "Failed.");

    return 0;
}
