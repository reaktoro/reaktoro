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
from numpy import *

# Construct the chemical system with its phases and species
editor = ChemicalEditor()
editor.addAqueousPhase('H2O NaCl CaCO3 MgCO3 CO2')
# editor.addAqueousPhase(['H2O(l)', 'H+', 'OH-', 'Na+', 'Cl-', 'Ca++', 'Mg++', 'HCO3-', 'CO2(aq)', 'CO3--'])
editor.addMineralPhase('Quartz')
editor.addMineralPhase('Calcite')
editor.addMineralPhase('Dolomite')

# Create the ChemicalSystem object using the configured editor
system = ChemicalSystem(editor)

print system
# 
# # Define the inicial condition of the reactive transport modeling problem
# problem_ic = EquilibriumProblem(system)
# problem_ic.add('H2O', 1.0, 'kg')
# problem_ic.add('NaCl', 0.7, 'mol')
# problem_ic.add('CaCO3', 10, 'mol')
# problem_ic.add('SiO2', 10, 'mol')
# 
# # Define the boundary condition of the reactive transport modeling problem
# problem_bc = EquilibriumProblem(system)
# problem_bc.add('H2O', 1.0, 'kg')
# problem_bc.add('NaCl', 0.90, 'mol')
# problem_bc.add('MgCl2', 0.05, 'mol')
# problem_bc.add('CaCl2', 0.01, 'mol')
# problem_bc.add('CO2', 0.75, 'mol')
# 
# # Calculate the equilibrium states for the initial and boundary conditions
# state_ic = equilibrate(problem_ic)
# state_bc = equilibrate(problem_bc)
# 
# # Scale the phases in the initial condition as required
# state_ic.scalePhaseVolume('Aqueous', 0.1, 'm3')
# state_ic.scalePhaseVolume('Quartz', 0.88, 'm3')
# state_ic.scalePhaseVolume('Calcite', 0.02, 'm3')
# 
# print state_bc
# print state_ic