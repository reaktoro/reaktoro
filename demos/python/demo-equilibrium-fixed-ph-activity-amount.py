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
editor.addAqueousPhase("H O Na Cl Ca Mg C")
editor.addGaseousPhase("H O C")

system = ChemicalSystem(editor)

problem = EquilibriumInverseProblem(system)
problem.add("H2O", 1, "kg")
problem.add("NaCl", 0.1, "mol")
problem.add("CaCl2", 2, "mmol")
problem.add("MgCl2", 4, "mmol")
problem.pH(3.0, "HCl")
problem.fixSpeciesAmount("CO2(g)", 1.0, "mol")
problem.fixSpeciesActivity("O2(g)", 0.20)

state = equilibrate(problem)

print state
