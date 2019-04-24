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

# Step 2: Add the phases of the chemical system
editor = ChemicalEditor()
editor.addAqueousPhase("H2O NaCl CaCO3 MgCO3")
#editor.addAqueousPhaseWithElementsOf("H2O NaCl CaCO3 MgCO3")
editor.addGaseousPhase(["H2O(g)", "CO2(g)"])
editor.addMineralPhase("Calcite")
editor.addMineralPhase("Magnesite")
editor.addMineralPhase("Dolomite")
editor.addMineralPhase("Halite")

# Step 3: Add the three reactions of the chemical system
editor.addMineralReaction("Calcite") \
    .setEquation("Calcite = Ca++ + CO3--") \
    .addMechanism("logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol") \
    .addMechanism("logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0") \
    .setSpecificSurfaceArea(10, "cm2/g")

editor.addMineralReaction("Magnesite") \
    .setEquation("Magnesite = Mg++ + CO3--") \
    .addMechanism("logk = -9.34 mol/(m2*s); Ea = 23.5 kJ/mol") \
    .addMechanism("logk = -6.38 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0") \
    .setSpecificSurfaceArea(10, "cm2/g")

editor.addMineralReaction("Dolomite") \
    .setEquation("Dolomite = Ca++ + Mg++ + 2*CO3--") \
    .addMechanism("logk = -7.53 mol/(m2*s); Ea = 52.2 kJ/mol") \
    .addMechanism("logk = -3.19 mol/(m2*s); Ea = 36.1 kJ/mol; a[H+] = 0.5") \
    .setSpecificSurfaceArea(10, "cm2/g")

# Step 4: Construct the chemical system
system = ChemicalSystem(editor)
reactions = ReactionSystem(editor)

# Step 5: Define the partition of the chemical system to the kinetics and equilibrium species
partition = Partition(system)
# Set the kinetics species
partition.setKineticSpecies(["Calcite", "Magnesite", "Dolomite"])

# Step 6: Define the chemical equilibrium problem
problem = EquilibriumProblem(system)
# Provide the partition of the equilibrium problem
problem.setPartition(partition)
# Add amounts of the species
problem.add("H2O", 1, "kg")
problem.add("NaCl", 1, "mol")
problem.add("CO2", 1, "mol")

# Step 7: Calculate the chemical equilibrium state
state0 = equilibrate(problem)

# Step 8: Setting the mass of mineral
state0.setSpeciesMass("Calcite", 100, "g")
state0.setSpeciesMass("Dolomite", 50, "g")

# Step 9: Model the kinetics' path of the reactions of the chemical system
path = KineticPath(reactions)
path.setPartition(partition)

# Step 10: Plot kinetics' paths
plot0 = path.plot()
plot0.x("time(units=hour)")
plot0.y("pH")
plot0.xlabel("Time [hour]")
plot0.ylabel("pH")
plot0.showlegend(False)

plot1 = path.plot()
plot1.x("time(units=hour)")
plot1.y("elementMolality(Ca)", "Ca")
plot1.y("elementMolality(Mg)", "Mg")
plot1.xlabel("Time [hour]")
plot1.ylabel("Concentration [molal]")
plot1.legend("right center")

plot2 = path.plot()
plot2.x("time(units=hour)")
plot2.y("phaseMass(Calcite units=grams)", "Calcite")
plot2.xlabel("Time [hour]")
plot2.ylabel("Mass [g]")

plot3 = path.plot()
plot3.x("time(units=hour)")
plot3.y("phaseMass(Dolomite units=grams)", "Dolomite")
plot3.xlabel("Time [hour]")
plot3.ylabel("Mass [g]")

# Step 11: Solve the kinetic paths
t0, t1 = 0.0, 25.0
path.solve(state0, t0, t1, "hours")
# Print the result state
print(state0)

