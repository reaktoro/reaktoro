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

# Step 2: Specify the phases in the chemical system and their species
editor = ChemicalEditor()
editor.addAqueousPhaseWithElementsOf("H2O NaCl CaCO3 MgCO3")
editor.addGaseousPhase(["H2O(g)", "CO2(g)"])
editor.addMineralPhase("Calcite")   # CaCO3
editor.addMineralPhase("Magnesite") # MgCO3
editor.addMineralPhase("Dolomite")  # CaMg(CO3)2
editor.addMineralPhase("Halite")    # HaCl

# Step 3: Define the kinetically-controlled reactions
editor.addMineralReaction("Calcite") \
    .setEquation("Calcite = Ca++ + CO3--") \
    .addMechanism("logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol") \
    .addMechanism("logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0") \
    .setSpecificSurfaceArea(10, "cm2/g")

editor.addMineralReaction("Dolomite") \
    .setEquation("Dolomite = Ca++ + Mg++ + 2*CO3--") \
    .addMechanism("logk = -7.53 mol/(m2*s); Ea = 52.2 kJ/mol") \
    .addMechanism("logk = -3.19 mol/(m2*s); Ea = 36.1 kJ/mol; a[H+] = 0.5") \
    .setSpecificSurfaceArea(10, "cm2/g")

editor.addMineralReaction("Magnesite") \
    .setEquation("Magnesite = Mg++ + CO3--") \
    .addMechanism("logk = -9.34 mol/(m2*s); Ea = 23.5 kJ/mol") \
    .addMechanism("logk = -6.38 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0") \
    .setSpecificSurfaceArea(10, "cm2/g")


# Step 4: Construct the chemical system
system = ChemicalSystem(editor)
reactions = ReactionSystem(editor)

# Step 5: Specify the equilibrium and kinetic species
partition = Partition(system)
partition.setKineticSpecies(["Calcite", "Magnesite", "Dolomite"])

# Step 6: Define the initial chemical equilibrium state
problem = EquilibriumProblem(partition)
problem.setTemperature(60, "celsius")
problem.setPressure(100, "bar")
problem.add("H2O", 1, "kg")
problem.add("NaCl", 1, "mol")
problem.add("CO2", 1, "mol")

# Step 7: Calculate the initial chemical equilibrium state
state0 = equilibrate(problem)
state0.output('demo-kineticpath-carbonates-co2-before-kinetics.txt')

# Step 8: Set the initial mass of the kinetic species
state0.setSpeciesMass("Calcite", 100, "g")
state0.setSpeciesMass("Dolomite", 50, "g")
state0.setSpeciesMass("Magnesite", 30, "g")

# Step 9: Create a kinetic path solver
path = KineticPath(reactions, partition)

# Step 10: Define properties to be output
output = path.output()
output.add("time(units=second)")
output.add("pH")
output.add("speciesMass(Calcite units=g)", "Calcite")
output.add("speciesMass(Dolomite units=g)", "Dolomite")
output.add("speciesMass(Magnesite units=g)", "Magnesite")
output.filename("demo-kineticpath-carbonates-co2-kinetics.txt")

# Step 11: Perform the kinetic path calculation
t0, t1 = 0.0, 25.0
path.solve(state0, t0, t1, "hours")

# Step 12: Output the final chemical state of the system
state0.output('demo-kineticpath-carbonates-co2-after-kinetics.txt')

# Step 13: Plotting the results of kinetic path calculation

import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
mpl.set_loglevel("critical")

# To load results from the outputfile, we use `loadtxt` function provided by the `numpy` package:

filearray = np.loadtxt("demo-kineticpath-carbonates-co2-kinetics.txt", skiprows=1)
data = filearray.T
[t_indx, ph_indx, calcite_indx, dolomite_indx, magnesite_indx] = np.arange(5)

plt.figure()
plt.plot(data[t_indx], data[calcite_indx], label="Calcite", color="C1")
plt.legend(loc='center right')
plt.xlabel("Time [s]")
plt.ylabel("Mass [g]")
plt.savefig("calcite-mass.png")

plt.figure()
plt.plot(data[t_indx], data[dolomite_indx], label="Dolomite", color="C2")
plt.legend(loc='center right')
plt.xlabel("Time [s]")
plt.ylabel("Mass [g]")
plt.savefig("dolomite-mass.png")

plt.figure()
plt.plot(data[t_indx], data[magnesite_indx], label="Magnesite", color="C3")
plt.legend(loc='center right')
plt.xlabel("Time [s]")
plt.ylabel("Mass [g]")
plt.savefig("magnesite-mass.png")

plt.figure()
plt.plot(data[t_indx], data[ph_indx], color="C4")
plt.xlabel("Time [s]")
plt.ylabel("pH [-]")
plt.savefig("ph.png")
