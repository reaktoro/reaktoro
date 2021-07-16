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

// -----------------------------------------------------------------------------
// 👏 Acknowledgements 👏
// -----------------------------------------------------------------------------
// This example was originally authored by:
//   • Allan Leal (16 July 2021)
//   • Svetlana Kyas (14 July 2021)
//
// and since revised by:
//   • Allan Leal (16 July 2021)
// -----------------------------------------------------------------------------

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    Database db("supcrtbl.yaml");

    AqueousPhase aqueousphase("H2O(aq) CO2(aq) CO3-2 Cl- H+ H2(aq) HCO3- Na+ NaCl(aq) NaOH(aq) O2(aq) OH- HCl(aq)");
    aqueousphase.setActivityModel(chain(
        ActivityModelHKF(),
        ActivityModelDrummond("CO2")
    ));

    GaseousPhase gaseousphase("CO2(g) H2O(g)");
    gaseousphase.setActivityModel(ActivityModelPengRobinson());

    Phases phases(db);
    phases.add(aqueousphase);
    phases.add(gaseousphase);

    ChemicalSystem system(phases);

    EquilibriumSpecs specs(system);
    specs.temperature();
    specs.pressure();
    specs.pH();

    EquilibriumOptions options;
    // options.optima.output.active = true;

    EquilibriumSolver solver(specs);
    solver.setOptions(options);

    EquilibriumConditions conditions(specs);
    conditions.temperature(60.0, "celsius");
    conditions.pressure(100.0, "bar");
    conditions.pH(4.0);
    conditions.startWith("H2O(aq)", 1.0, "kg");
    conditions.startWith("Na+",     1.0, "mol");
    conditions.startWith("Cl-",     1.0, "mol");
    conditions.startWith("CO2(g)", 10.0, "mol");

    ChemicalState state(system);

    solver.solve(state, conditions);

    std::cout << state << std::endl;

    ChemicalProps props(state);
    std::cout << props << std::endl;

    AqueousProps aprops(state);

    std::cout << aprops << std::endl;

    return 0;
}
