// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021
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
    Database db("supcrt98.yaml");

    // Create an aqueous phase without organic species
    AqueousPhase aqphase(speciate("H O C Na Cl"), exclude("organic"));
    aqphase.setActivityModel(chain(
            ActivityModelHKF(),
            ActivityModelDrummond("CO2")
    ));

    // Create a gaseous phase
    GaseousPhase gasphase("CO2(g)");
    gasphase.setActivityModel(ActivityModelPengRobinson());

    std::cout << "Phases with organics" << std::endl;

    // Collecting all above-defined phases
    Phases phases(db);
    phases.add(aqphase);
    phases.add(gasphase);

    // Construct the chemical system
    ChemicalSystem system(phases);
    // Define initial equilibrium state
    ChemicalState state(system);
    state.setTemperature(25.0, "celsius");
    state.setPressure(1.0, "bar");
    state.setSpeciesMass("H2O(aq)", 1.0, "kg");
    state.setSpeciesAmount("CO2(g)", 10.0, "mol");
    state.setSpeciesAmount("Na+", 4.0, "mol");
    state.setSpeciesAmount("Cl-", 4.0, "mol");

    // Define equilibrium solver and equilibrate given initial state
    EquilibriumSolver solver(system);
    solver.solve(state);

    // Obtain species composition from the equilibrated state
    const auto n = state.speciesAmounts();

    // Print the species and theirs amounts
    std::cout << std::setw(20) << "Species"
              << std::setw(20) << "Amount" << std::endl;

    for(auto i = 0; i < n.size(); ++i)
    {
        // Print only species with nonzero amounts
        if (n[i] > 1e-16)
            std::cout << std::setw(20) << system.species(i).name()
                      << std::setw(20) << n[i] << std::endl;
    }

    return 0;
}