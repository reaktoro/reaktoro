// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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
    aqueousphase.setActivityModel(chain(
        ActivityModelHKF(),
        ActivityModelDrummond("CO2")
    ));

    // Define gaseous phase
    GaseousPhase gaseousphase("CO2(g)");
    gaseousphase.setActivityModel(ActivityModelPengRobinson());

    // Initialize phases with aqueous and gaseous phase
    Phases phases(db);
    phases.add(aqueousphase);
    phases.add(gaseousphase);

    // Construct the chemical system
    ChemicalSystem system(phases);

    // Define initial equilibrium state
    ChemicalState state(system);
    state.setTemperature(25.0, "celsius");
    state.setPressure(1.0, "bar");
    state.setSpeciesMass("H2O"     , 1.00, "kg");
    state.setSpeciesAmount("CO2(g)", 10.0, "mol");
    state.setSpeciesAmount("Na+"   , 4.00, "mol");
    state.setSpeciesAmount("Cl-"   , 4.00, "mol");

    // Define equilibrium solver and equilibrate given initial state
    EquilibriumSolver solver(system);
    solver.solve(state);

    // Output chemical state to the txt-file
    state.output("state.txt");

    return 0;
}
