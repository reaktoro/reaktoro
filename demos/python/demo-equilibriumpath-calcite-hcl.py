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
editor.addAqueousPhase("H O Ca C Cl")
editor.addMineralPhase("Calcite")

system = ChemicalSystem(editor)

problem1 = EquilibriumProblem(system)
problem1.setTemperature(30.0, "celsius")
problem1.setPressure(1.0, "bar")
problem1.add("H2O", 1, "kg")
problem1.add("CaCO3", 100, "g")

problem2 = EquilibriumProblem(system)
problem2.setTemperature(30.0, "celsius")
problem2.setPressure(1.0, "bar")
problem2.add("H2O", 1, "kg")
problem2.add("CaCO3", 100, "g")
problem2.add("HCl", 1, "mmol")

state1 = equilibrate(problem1)
state2 = equilibrate(problem2)

path = EquilibriumPath(system)

plot1 = path.plot()
plot1.x("elementAmount(Cl units=mmol)")
plot1.y("pH", "pH")
plot1.xlabel("HCl [mmol]")
plot1.ylabel("pH")
plot1.showlegend(False)

plot2 = path.plot()
plot2.x("pH")
plot2.y("Ca", "elementMolality(Ca units=mmolal)")
plot2.xlabel("pH")
plot2.ylabel("Concentration [mmolal]")

plot3 = path.plot()
plot3.x("pH")
plot3.y("CO2(aq)", "speciesMolality(CO2(aq) units=mmolal)")
plot3.y("CO3--",   "speciesMolality(CO3-- units=mmolal)")
plot3.xlabel("pH")
plot3.ylabel("Concentration [mmolal]")
plot3.legend("right bottom")

output = path.output()
output.file("result.txt")
output.add("Cl [mmol]", "elementAmount(Cl units=mmol)")
output.add("Ca [molal]", "elementMolality(Ca) pH")
output.add("pH")

path.solve(state1, state2)
