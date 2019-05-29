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

from reaktoro import *

editor = ChemicalEditor()
editor.addAqueousPhaseWithElementsOf("H2O NaCl CaCO3")
editor.addGaseousPhase(["H2O(g)", "CO2(g)"])
editor.addMineralPhase("Calcite")

system = ChemicalSystem(editor)

problem = EquilibriumInverseProblem(system)
problem.add("H2O", 1, "kg")
problem.add("NaCl", 0.1, "mol")
problem.fixPhaseVolume("Gaseous", 0.2, "m3", "CO2")
problem.fixPhaseVolume("Aqueous", 0.3, "m3", "1 kg H2O; 0.1 mol NaCl")
problem.fixPhaseVolume("Calcite", 0.5, "m3", "CaCO3")

state = equilibrate(problem)

print(state)
