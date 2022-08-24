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

// -----------------------------------------------------------------------------
// 👏 Acknowledgements 👏
// -----------------------------------------------------------------------------
// This example was originally authored by:
//   • Allan Leal (4 August 2021)
// -----------------------------------------------------------------------------

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    SupcrtDatabase db("supcrtbl");

    AqueousPhase solution(speciate("Na Cl C"), exclude("organic"));
    solution.setActivityModel(chain(
        ActivityModelHKF(),
        ActivityModelDrummond("CO2")
    ));

    GaseousPhase gases("CO2(g) H2O(g)");
    gases.setActivityModel(ActivityModelPengRobinson());

    ChemicalSystem system(db, solution, gases);

    EquilibriumSpecs specs(system);
    specs.temperature();
    specs.phaseAmount("GaseousPhase");

    EquilibriumSolver solver(specs);

    ChemicalState state(system);
    state.setTemperature(50.0, "celsius");
    state.setPressure(300.0, "bar");
    state.set("H2O(aq)", 1.0, "kg");
    state.set("Na+",     1.0, "mol");
    state.set("Cl-",     1.0, "mol");
    state.set("CO2(aq)", 1.0, "mol");

    EquilibriumConditions conditions(specs);
    conditions.temperature(50.0, "celsius");
    conditions.phaseAmount("GaseousPhase", 1e-10, "mol");
    conditions.setLowerBoundPressure(1.0, "bar");
    conditions.setUpperBoundPressure(1000.0, "bar");

    solver.solve(state, conditions);

    std::cout << "Pressure: " << state.pressure() * 1e-5 << " bar" << std::endl;

    return 0;
}
