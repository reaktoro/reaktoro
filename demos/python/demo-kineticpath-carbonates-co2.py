# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright (C) 2014-2015 Allan Leal
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

from reaktoro import *

editor = ChemicalEditor()
editor.addAqueousPhase("H2O NaCl CaCO3 MgCO3")
editor.addGaseousPhase(["H2O(g)", "CO2(g)"])
editor.addMineralPhase("Calcite")
editor.addMineralPhase("Magnesite")
editor.addMineralPhase("Dolomite")
editor.addMineralPhase("Halite")

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

system = ChemicalSystem(editor)
reactions = ReactionSystem(editor)

partition = Partition(system)
partition.setKineticSpecies(["Calcite", "Magnesite", "Dolomite"])

problem = EquilibriumProblem(system)
problem.setPartition(partition)
problem.add("H2O", 1, "kg")
problem.add("NaCl", 1, "mol")
problem.add("CO2", 1, "mol")

state0 = equilibrate(problem)

state0.setSpeciesMass("Calcite", 100, "g")
state0.setSpeciesMass("Dolomite", 50, "g")

path = KineticPath(reactions)
path.setPartition(partition)

plot0 = path.plot()
plot0.x("time(units=hour)")
plot0.y("pH", "pH")
plot0.xlabel("Time [hour]")
plot0.ylabel("pH")
plot0.showlegend(False)

plot1 = path.plot()
plot1.x("time(units=hour)")
plot1.y("Ca", "elementMolality(Ca)")
plot1.y("Mg", "elementMolality(Mg)")
plot1.xlabel("Time [hour]")
plot1.ylabel("Concentration [molal]")
plot1.legend("right center")

plot2 = path.plot()
plot2.x("time(units=hour)")
plot2.y("Calcite", "phaseMass(Calcite units=grams)")
plot2.xlabel("Time [hour]")
plot2.ylabel("Mass [g]")

plot3 = path.plot()
plot3.x("time(units=hour)")
plot3.y("Dolomite", "phaseMass(Dolomite units=grams)")
plot3.xlabel("Time [hour]")
plot3.ylabel("Mass [g]")

path.solve(state0, 0, 25, "hours")

print state0

