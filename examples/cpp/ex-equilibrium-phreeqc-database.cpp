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
//   • Svetlana Kyas (14 July 2021)
//
// and since revised by:
//   • Allan Leal (28 August 2023)
//     - Using ActivityModelPhreeqc instead of ActivityModelHKF for aqueous phase.
//     - Using ActivityModelPengRobinsonPhreeqc instead of ActivityModelPengRobinson for gaseous phase
//   • Allan Leal (16 July 2021)
// -----------------------------------------------------------------------------

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    // Initialize a thermodynamic database
    PhreeqcDatabase db("phreeqc.dat");

    // Define aqueous phase
    AqueousPhase aqueousphase(speciate("H O C Na Cl"));
    aqueousphase.set(ActivityModelPhreeqc(db));

    // Define gaseous phase
    GaseousPhase gaseousphase("CO2(g)");
    gaseousphase.set(ActivityModelPengRobinsonPhreeqc());

    // Initialize phases with aqueous and gaseous phase
    Phases phases(db);
    phases.add(aqueousphase);
    phases.add(gaseousphase);

    // Construct the chemical system
    ChemicalSystem system(phases);

    // Define initial equilibrium state
    ChemicalState state(system);
    state.temperature(25.0, "celsius");
    state.pressure(1.0, "bar");
    state.set("H2O"     , 1.00, "kg");
    state.set("CO2(g)", 10.0, "mol");
    state.set("Na+"   , 4.00, "mol");
    state.set("Cl-"   , 4.00, "mol");

    // Create an equilibrium solver and equilibrate given initial state
    EquilibriumSolver solver(system);
    solver.solve(state);

    // Output the chemical state to a text file
    state.output("state.txt");

    return 0;
}
