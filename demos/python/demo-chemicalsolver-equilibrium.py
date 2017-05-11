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


# FIXME



# from reaktoro import *
# from numpy import *

# npoints = 10

# editor = ChemicalEditor()
# editor.addAqueousPhase("H O Na Cl C Ca Mg")
# editor.addGaseousPhase(["H2O(g)", "CO2(g)"])
# editor.addMineralPhase("Calcite")
# editor.addMineralPhase("Dolomite")

# system = ChemicalSystem(editor)

# composition = EquilibriumCompositionProblem(system)
# composition.setAqueousComposition("1 molal NaCl")
# composition.setGaseousComposition("CO2")
# composition.setSolidComposition("0.1 Calcite; 0.9 Dolomite")
# composition.setAqueousSaturation(0.8)
# composition.setGaseousSaturation(0.2)
# composition.setPorosity(0.3)

# state = equilibrate(composition)

# solver = ChemicalSolver(system, npoints)

# Ee = system.numElements()

# T = array([state.temperature()] * npoints)
# P = array([state.pressure()] * npoints)
# be = zeros((Ee, npoints))

# for i in xrange(npoints):
#     be[:, i] = state.elementAmounts().array()

# print T
# print P
# print be

# print shape(T)
# print shape(P)
# print shape(be)

# solver.equilibrate(T, P, be)

# for state in solver.states():
#     print state

# print "porosity = \n", solver.porosity()
# print "densities[0] = \n", solver.fluidDensities()[0]
# print "densities[1] = \n", solver.fluidDensities()[1]
# print "saturations[0] = \n", solver.fluidSaturations()[0]
# print "saturations[1] = \n", solver.fluidSaturations()[1]
