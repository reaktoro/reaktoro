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
editor.addAqueousPhase("H O Na Cl Ca C")
editor.addMineralPhase("Calcite")

system = ChemicalSystem(editor)

problem = EquilibriumInverseProblem(system)
problem.add("H2O", 1, "kg")
problem.add("NaCl", 0.1, "mol")
problem.pH(8.0, "HCl", "NaOH")
problem.fixSpeciesAmount("Calcite", 1, "mol")

state = equilibrate(problem)

print state
