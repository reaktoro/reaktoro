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

db = Database('resources/supcrt98-william-smith-equilibrium-challenge.xml')

editor = ChemicalEditor(db)
phase = editor.addGaseousPhase([
    'CO(g)',
    'CO2(g)',
    'Cl2(g)',
    'HCl(g)',
    'H2(g)',
    'H2O(g)',
    'O2(g)']).setChemicalModelIdeal()

system = ChemicalSystem(editor)

problem = EquilibriumProblem(system)
problem.setTemperature(300, "K")
problem.setPressure(1, "bar")
problem.add("CO",  1.0, "mol")
problem.add("Cl2", 0.5, "mol")
problem.add("H2",  1.5, "mol")
problem.add("O2",  1.0, "mol")

state = equilibrate(problem)

state.output('demo-equilibrium-william-smith-challenge.txt')
