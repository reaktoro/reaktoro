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

# Step 2: Define a chemical system with an aqueous and a mineral phase
editor = ChemicalEditor()
editor.addAqueousPhaseWithElementsOf("H2O HCl CaCO3")
editor.addMineralPhase("Calcite")

# Step 3: Define mineral reaction for Calcite
editor.addMineralReaction("Calcite") \
    .setEquation("Calcite = Ca++ + CO3--") \
    .addMechanism("logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol") \
    .addMechanism("logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0") \
    .setSpecificSurfaceArea(10, "cm2/g")

# Step 4: Create the ChemicalSystem and ReactionSystem instances
system = ChemicalSystem(editor)
reactions = ReactionSystem(editor)

# Step 5: Specify the equilibrium and kinetic species
partition = Partition(system)
partition.setKineticSpecies(["Calcite"])

# Step 6: Define the chemical equilibrium problem for the equilibrium partition
problem = EquilibriumProblem(partition)
problem.setTemperature(30, "celsius")
problem.setPressure(1, "bar")
problem.add("H2O", 1, "kg")
problem.add("HCl", 1, "mmol")

# Step 7: Calculate the chemical equilibrium state in the equilibrium partition
state0 = equilibrate(problem)
state0.output('demo-kineticpath-calcite-hcl-before-kinetics.txt')

# Step 8: Set the initial mass of the kinetic species
state0.setSpeciesMass("Calcite", 100, "g")

# Step 9: Define the kinetic path problem
path = KineticPath(reactions, partition)

# Step 10: Define properties to be output
output = path.output()
output.filename("result.txt")
output.add("time(units=minute)")
output.add("elementMolality(Ca units=mmolal)", "Ca")
output.add("phaseMass(Calcite units=g)", "Calcite")
output.add("speciesMolality(Ca++ units=mmolal)", "Ca++")
output.add("speciesMolality(HCO3- units=mmolal)", "HCO3-")
output.add("pH")

# Step 10: Perform the kinetic path calculation
t0, t1 = 0.0, 5.0
path.solve(state0, t0, t1, "minute")

# Step 11: Output the final state of the chemical system
state0.output('demo-kineticpath-calcite-hcl-after-kinetics')

# Step 12: Plot different properties of the chemical system during kinetics

import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
mpl.set_loglevel("critical")

filearray = np.loadtxt("result.txt", skiprows=1)
data = filearray.T
[t_indx, ca_indx, calcite_indx, ca2_indx, hco3_indx, ph_indx] = np.arange(6)

plt.figure()
plt.plot(data[t_indx], data[ca_indx], label="Ca", color="C1")
plt.xlabel("Time [minute]")
plt.ylabel("Concentration [mmolal]")
plt.legend(loc='center right')
plt.savefig("ca-concentration.png")

plt.figure()
plt.plot(data[t_indx], data[calcite_indx], label="Calcite", color="C2")
plt.xlabel("Time [minute]")
plt.ylabel("Phase mass [g]")
plt.legend(loc='center right')
plt.savefig("calcite-phase-mass.png")

plt.figure()
plt.plot(data[t_indx], data[ph_indx], color="C3")
plt.xlabel("Time [minute]")
plt.ylabel("pH [-]")
plt.savefig("ph.png")

plt.figure()
plt.plot(data[t_indx], data[ca2_indx], label="Ca++", color="C4")
plt.xlabel("Time [minute]")
plt.ylabel("Concentration [mmolal]")
plt.legend(loc='center right')
plt.savefig("aqueous-species-concentration.png")
