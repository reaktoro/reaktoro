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

# Auxiliary time related constants
second = 1
minute = 60
hour = 60 * minute
day = 24 * hour
year = 365 * day

# Parameters for the reactive transport simulation
nsteps = 100      # the number of steps in the reactive transport simulation
ncells = 100      # the number of cells in the discretization
xl = 0.0          # the x-coordinate of the left boundary
xr = 100.0        # the x-coordinate of the right boundary
D  = 1.0e-9       # the diffusion coefficinet (in units of m2/s)
v  = 1.0/day      # the fluid pore velocity (in units of m/s)
dt = 0.5*day      # the time step (in units of s)
T = 60.0          # the temperature (in units of degC)
P = 100           # the pressure (in units of bar)

# Construct the chemical system with its phases and species
editor = ChemicalEditor()
editor.addAqueousPhase('H2O NaCl CaCO3 MgCO3 CO2')
editor.addAqueousPhase(['H2O(l)', 'H+', 'OH-', 'Na+', 'Cl-', 'Ca++', 'Mg++', 'HCO3-', 'CO2(aq)', 'CO3--'])
editor.addMineralPhase('Quartz')
editor.addMineralPhase('Calcite')
editor.addMineralPhase('Dolomite')

# Create the ChemicalSystem object using the configured editor
system = ChemicalSystem(editor)

# Define the inicial condition of the reactive transport modeling problem
problem_ic = EquilibriumProblem(system)
problem_ic.setTemperature(T, 'celsius')
problem_ic.setPressure(P, 'bar')
problem_ic.add('H2O', 1.0, 'kg')
problem_ic.add('NaCl', 0.7, 'mol')
problem_ic.add('CaCO3', 10, 'mol')
problem_ic.add('SiO2', 10, 'mol')
 
# Define the boundary condition of the reactive transport modeling problem
problem_bc = EquilibriumProblem(system)
problem_bc.setTemperature(T, 'celsius') 
problem_bc.setPressure(P, 'bar')
problem_bc.add('H2O', 1.0, 'kg')
problem_bc.add('NaCl', 0.90, 'mol')
problem_bc.add('MgCl2', 0.05, 'mol')
problem_bc.add('CaCl2', 0.01, 'mol')
problem_bc.add('CO2', 0.75, 'mol')
 
# Calculate the equilibrium states for the initial and boundary conditions
state_ic = equilibrate(problem_ic)
state_bc = equilibrate(problem_bc)
 
# Scale the phases in the initial condition as required
state_ic.scalePhaseVolume('Aqueous', 0.1, 'm3')
state_ic.scalePhaseVolume('Quartz', 0.88, 'm3')
state_ic.scalePhaseVolume('Calcite', 0.02, 'm3')

# Create the mesh for the column
mesh = Mesh(ncells, xl, xr)

# Create a chemical field object with every cell having state given by state_ic
field = ChemicalField(mesh.numCells(), state_ic)

# Define the reactive transport modeling
rt = ReactiveTransportSolver(system)
rt.setMesh(mesh)
rt.setVelocity(v)
rt.setDiffusionCoeff(D)
rt.setBoundaryState(state_bc)
rt.setTimeStep(dt)

# Define the quantities that should be output for every cell, every time step
output = rt.output()
output.filename('reativetransport.txt')
output.add('pH')
output.add('speciesMolality(H+)')
output.add('speciesMolality(Ca++)')
output.add('speciesMolality(Mg++)')
output.add('speciesMolality(HCO3-)')
output.add('speciesMolality(CO2(aq))')
output.add('phaseVolume(Calcite)')
output.add('phaseVolume(Dolomite)')

rt.initialize(field)

t = 0.0
step = 0

while step <= nsteps:
    # Print the progress of the simulation
    print "Progress: {}/{} steps, {} min".format(step, nsteps, t/minute)
    
    # Perform one reactive transport time step
    rt.step(field)
    
    # Increment time step and number of time steps
    t += dt
    step += 1


# reactivetransport = ReactiveTransportSolver(system)
# reactivetransport.setMesh(mes)
