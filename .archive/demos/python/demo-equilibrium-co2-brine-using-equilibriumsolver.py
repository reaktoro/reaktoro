# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright (C) 2014-2020 Allan Leal
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library. If not, see <http://www.gnu.org/licenses/>.

from reaktoro import *

# Create an object of Database class to use in the initialization of the
# chemical system.
database = Database("supcrt98.xml")

# Define which phases and species the chemical system should have using a
# ChemicalEditor object.
editor = ChemicalEditor(database)
editor.addAqueousPhaseWithElementsOf("H2O NaCl CO2")
editor.addGaseousPhase(["H2O(g)", "CO2(g)"])
editor.addMineralPhase("Halite")

# Initialize the chemical system using the definition provided to the
# chemical editor.
system = ChemicalSystem(editor)

# Define an equilibrium problem with given temperature, pressure, and
# amounts of compounds.
problem = EquilibriumProblem(system)
problem.setTemperature(60, "celsius")
problem.setPressure(300, "bar")
problem.add("H2O", 1, "kg")
problem.add("CO2", 100, "g")
problem.add("NaCl", 0.1, "mol")

# Get the temperature, pressure, and mole amounts of the elements.
T = problem.temperature()
P = problem.pressure()
b = problem.elementAmounts()

# Create an object of EquilibriumSolver class that can be reused many times.
solver = EquilibriumSolver(system)

# Create an object of ChemicalState class to store the equilibrium
# state of the system.
state = ChemicalState(system)

# Solve the equilibrium state with given (T, P, b) inputs.
solver.solve(state, T, P, b)

# Print the calculated chemical equilibrium state.
print(state)

# Calculate the new equilibrium state when temperature is increased.
# Use the previous equilibrium state as an initial guess for improved
# performance.
solver.solve(state, T + 10.0, P, b)

# Print the new calculated chemical equilibrium state.
print(state)
