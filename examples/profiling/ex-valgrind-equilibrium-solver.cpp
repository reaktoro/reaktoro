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
    SupcrtDatabase db("supcrtbl");

    AqueousPhase solution("H2O(aq) H+ OH- Na+ Cl- HCO3- CO3-2 CO2(aq)");
    solution.setActivityModel(chain(
        ActivityModelHKF(),
        ActivityModelDrummond("CO2")
    ));

    GaseousPhase gases("CO2(g) H2O(g)");
    gases.setActivityModel(ActivityModelPengRobinson());

    ChemicalSystem system(db, solution, gases);

    EquilibriumOptions options;
    options.optima.newtonstep.linearsolver.method = Optima::LinearSolverMethod::Rangespace;
    options.hessian = GibbsHessian::ApproxDiagonal;

    EquilibriumSolver solver(system);
    solver.setOptions(options);

    ChemicalState state(system);
    state.temperature(60.0, "celsius");
    state.pressure(100.0, "bar");
    state.set("H2O(aq)", 1.0, "kg");
    state.set("Na+",     1.0, "mol");
    state.set("Cl-",     1.0, "mol");
    state.set("CO2(g)", 10.0, "mol");

    CALLGRIND_TOGGLE_COLLECT; // turn on
    auto result = solver.solve(state);
    CALLGRIND_TOGGLE_COLLECT; // turn off

    errorif(result.failed(), "Equilibrium calculation failed.");

    return 0;
}
