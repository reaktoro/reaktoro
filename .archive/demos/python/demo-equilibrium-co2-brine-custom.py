# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright (C) 2014-2021 Allan Leal
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

editor = ChemicalEditor()
editor.addAqueousPhase(["H2O(l)", "H+", "OH-", "Na+", "Cl-", "HCO3-", "CO2(aq)", "CO3--"]) \
    .setActivityModelDrummondCO2()
editor.addGaseousPhase(["H2O(g)", "CO2(g)"]). \
    setChemicalModelSpycherPruessEnnis()
editor.addMineralPhase("Halite")

system = ChemicalSystem(editor)

problem = EquilibriumProblem(system)
problem.setTemperature(60, "celsius")
problem.setPressure(300, "bar")
problem.add("H2O", 1, "kg")
problem.add("CO2", 100, "g")
problem.add("NaCl", 1, "mol")

state = equilibrate(problem)

print(state)
