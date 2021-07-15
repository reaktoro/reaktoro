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

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    // Initialize a thermodynamic database
    SupcrtDatabase db("supcrt98.xml");

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
    state.setSpeciesMass("H2O(l)", 1.0, "kg");
    state.setSpeciesAmount("CO2(g)", 10.0, "mol");
    state.setSpeciesAmount("Na+", 4.0, "mol");
    state.setSpeciesAmount("Cl-", 4.0, "mol");

    // Define equilibrium solver and equilibrate given initial state
    EquilibriumSolver solver(system);
    solver.solve(state);

    // Output temperature, pressure, and species amounts of the chemical state
    std::cout << state << std::endl;

    return 0;
}
