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
editor.addAqueousPhase("H2O HCl CaCO3")
#editor.addAqueousPhaseWithElementsOf("H2O HCl CaCO3")
editor.addMineralPhase("Calcite")

# Step 3: Define a mineral reactions
editor.addMineralReaction("Calcite") \
    .setEquation("Calcite = Ca++ + CO3--") \
    .addMechanism("logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol") \
    .addMechanism("logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0") \
    .setSpecificSurfaceArea(10, "cm2/g")

# Step 4: Create the ChemicalSystem and ReactionSystem instances
system = ChemicalSystem(editor)
reactions = ReactionSystem(editor)

# Step 5: Define the partition of the chemical system to the kinetics and equilibrium species
partition = Partition(system)
# Set the kinetics species
partition.setKineticPhases(["Calcite"])

# Step 6: Define the chemical equilibrium problem
problem = EquilibriumProblem(system)
# Provide the partition of the equilibrium problem
problem.setPartition(partition)
problem.add("H2O", 1, "kg")
problem.add("HCl", 1, "mmol")

# Step 7: Calculate the chemical equilibrium state
state0 = equilibrate(problem)
print(state0)

# Step 8: Specify the mass of the species
state0.setSpeciesMass("Calcite", 100, "g")

# Step 9: Model the kinetics' path of the reactions of the chemical system
path = KineticPath(reactions)
path.setPartition(partition)

# Step 10: Plotting of evolution of different properties of the chemical system
plot1 = path.plot()
plot1.x("time(units=minute)")
plot1.y("elementMolality(Ca units=mmolal)", "Ca")
plot1.xlabel("Time [minute]")
plot1.ylabel("Concentration [mmolal]")
plot1.legend("right center")

plot2 = path.plot()
plot2.x("time(units=minute)")
plot2.y("phaseMass(Calcite units=g)", "Calcite")
plot2.xlabel("Time [minute]")
plot2.ylabel("Mass [g]")

plot3 = path.plot()
plot3.x("time(units=minute)")
plot3.y("pH")
plot3.xlabel("Time [minute]")
plot3.ylabel("pH")
plot3.legend("right center")

plot4 = path.plot()
plot4.x("time(units=minute)")
plot4.y("speciesMolality(Ca++ units=mmolal)", "Ca++")
plot4.y("speciesMolality(HCO3- units=mmolal)", "HCO3-")
plot4.xlabel("Time [minute]")
plot4.ylabel("Concentration [mmolal]")
plot4.legend("right center")

# Step 11: Solve the kinetic paths
t0, t1 = 0.0, 5.0
path.solve(state0, t0, t1, "minute")
print(state0)
