# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright (C) 2014-2018 Allan Leal
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

# Step 1: Import the reaktoro Python package
from reaktoro import *

# Step 2: Initialize a thermodynamic database
db = Database("supcrt98.xml")

# Step 3: Define the chemical system
editor = ChemicalEditor(db)
editor.addAqueousPhaseWithElements("H O Na Cl C")
editor.addGaseousPhase(["CO2(g)"])

# Step 4: Construct the chemical system
system = ChemicalSystem(editor)

# Step 5: Define the chemical equilibrium problem
problem = EquilibriumProblem(system)
problem.setTemperature(60, "celsius")
problem.setPressure(100, "bar")
problem.add("H2O", 1.0, "kg")
problem.add("NaCl", 1.0, "mol")
problem.add("CO2", 10.0, "mol")

# Step 6: Calculate the chemical equilibrium state
state = equilibrate(problem)

# Step 7: Output the calculated chemical state to a file
state.output("result.txt")

# Step 8: Print the amounts of some aqueous species
print("Amount of CO2(aq):", state.speciesAmount("CO2(aq)"))
print("Amount of HCO3-:", state.speciesAmount("HCO3-"))
print("Amount of CO3--:", state.speciesAmount("CO3--"))

# Step 9: Print the amounts of element C in both aqueous and gaseous phases
print("Amount of C in aqueous phase:", state.elementAmountInPhase("C", "Aqueous"))
print("Amount of C in gaseous phase:", state.elementAmountInPhase("C", "Gaseous"))
