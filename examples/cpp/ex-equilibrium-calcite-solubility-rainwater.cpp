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

    // Create an aqueous phase automatically selecting all species with provided elements
    AqueousPhase aqueousphase(speciate("H O C Ca Mg K Cl Na S N"), exclude("organic"));
    aqueousphase.setActivityModel(chain(
        ActivityModelHKF(),
        ActivityModelDrummond("CO2")
    ));

    // Create a gaseous phase
    GaseousPhase gaseousphase("CO2(g)");
    gaseousphase.setActivityModel(ActivityModelPengRobinson());

    // Create a mineral phase
    MineralPhase mineralphase("Calcite");

    // Collecting all above-defined phases
    Phases phases(db);
    phases.add(aqueousphase);
    phases.add(gaseousphase);
    phases.add(mineralphase);

    // Construct the chemical system
    ChemicalSystem system(phases);

    // Set temperature and pressure
    double T = 25.0;
    double P = 1.0;
    double m0Calcite = 10.0;

    // Define initial equilibrium state
    ChemicalState state(system);
    state.setTemperature(T, "celsius");
    state.setPressure(P, "bar");

    // Specify rainwater composition
    state.setSpeciesMass("H2O(aq)", 1.0, "kg");
    state.setSpeciesMass("Na+", 2.05, "mg");
    state.setSpeciesMass("K+", 0.35, "mg");
    state.setSpeciesMass("Ca+2", 1.42, "mg");
    state.setSpeciesMass("Mg+2", 0.39, "mg");
    state.setSpeciesMass("Cl-", 3.47, "mg");
    state.setSpeciesMass("S2O4-2", 2.19, "mg");
    state.setSpeciesMass("NO3-", 0.27, "mg");
    state.setSpeciesMass("NH4+", 0.41, "mg");
    // Specify the amount of calcite
    state.setSpeciesAmount("Calcite", m0Calcite, "mol");

    // Define equilibrium solver and equilibrate given initial state
    EquilibriumSolver solver(system);
    solver.solve(state);

    // Output temperature, pressure, and species amounts of the chemical state
    std::cout << "Solubility of the Calcite in the rainwater : " << (m0Calcite - state.speciesAmount("Calcite")) / state.speciesMass("H2O(aq)") << std::endl;

    return 0;
}