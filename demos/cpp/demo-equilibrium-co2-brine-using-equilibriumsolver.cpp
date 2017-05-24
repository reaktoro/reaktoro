// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    // Create an object of Database class to use in the initialization of the chemical system.
    Database database("supcrt98.xml");

    // Define which phases and species the chemical system should have using a ChemicalEditor object.
    ChemicalEditor editor(database);
    editor.addAqueousPhase("H2O NaCl CO2");
    editor.addGaseousPhase({"H2O(g)", "CO2(g)"});
    editor.addMineralPhase("Halite");

    // Initialize the chemical system using the definition provided to the chemical editor.
    ChemicalSystem system(editor);

    // Define an equilibrium problem with given temperature, pressure, and amounts of compounds.
    EquilibriumProblem problem(system);
    problem.setTemperature(60, "celsius");
    problem.setPressure(300, "bar");
    problem.add("H2O", 1, "kg");
    problem.add("CO2", 100, "g");
    problem.add("NaCl", 0.1, "mol");

    // Get the temperature, pressure, and mole amounts of the elements.
    const double T = problem.temperature();
    const double P = problem.pressure();
    const Vector b = problem.elementAmounts();

    // Create an object of EquilibriumSolver class that can be reused many times.
    EquilibriumSolver solver(system);

    // Create an object of ChemicalState class to store the equilibrium state of the system.
    ChemicalState state(system);

    // Solve the equilibrium state with given (T, P, b) inputs.
    solver.solve(state, T, P, b);

    // Print the calculated chemical equilibrium state.
    std::cout << state << std::endl;

    // Calculate the new equilibrium state when temperature is increased.
    // Use the previous equilibrium state as an initial guess for improved performance.
    solver.solve(state, T + 10.0, P, b);

    // Print the new calculated chemical equilibrium state.
    std::cout << state << std::endl;
}
