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

# Step 1: Import the reaktoro Python package
from reaktoro import *
import os

# Step 2: Initialise auxiliary time-related constants
second = 1
minute = 60
hour = 60 * minute
day = 24 * hour
year = 365 * day

# Step 3: Define parameters for the reactive transport simulation
xl = 0.0          # the x-coordinate of the left boundary
xr = 100.0        # the x-coordinate of the right boundary
nsteps = 100      # the number of steps in the reactive transport simulation
ncells = 100      # the number of cells in the discretization
D  = 1.0e-9       # the diffusion coefficient (in units of m2/s)
v  = 1.0/day      # the fluid pore velocity (in units of m/s)
dt = 0.5*day      # the time step (in units of s)
T = 60.0          # the temperature (in units of degC)
P = 100           # the pressure (in units of bar)

# Step 4: Construct the chemical system with its phases and species
editor = ChemicalEditor()
editor.addAqueousPhaseWithElementsOf('H2O NaCl CaCl2 MgCl2 CO2')
editor.addMineralPhase('Quartz')
editor.addMineralPhase('Calcite')
editor.addMineralPhase('Dolomite')

# Step 5: Create the ChemicalSystem object using the configured editor
system = ChemicalSystem(editor)

# Step 6: Define the initial condition of the reactive transport modeling problem
problem_ic = EquilibriumProblem(system)
problem_ic.setTemperature(T, 'celsius')
problem_ic.setPressure(P, 'bar')
problem_ic.add('H2O', 1.0, 'kg')
problem_ic.add('NaCl', 0.7, 'mol')
problem_ic.add('CaCO3', 10, 'mol')
problem_ic.add('SiO2', 10, 'mol')
 
# Step 7: Define the boundary condition of the reactive transport modeling problem
problem_bc = EquilibriumProblem(system)
problem_bc.setTemperature(T, 'celsius') 
problem_bc.setPressure(P, 'bar')
problem_bc.add('H2O', 1.0, 'kg')
problem_bc.add('NaCl', 0.90, 'mol')
problem_bc.add('MgCl2', 0.05, 'mol')
problem_bc.add('CaCl2', 0.01, 'mol')
problem_bc.add('CO2', 0.75, 'mol')
 
# Step 8: Calculate the equilibrium states for the initial and boundary conditions
state_ic = equilibrate(problem_ic)
state_bc = equilibrate(problem_bc)
 
# Step 9: Scale the volumes of the phases in the initial condition
state_ic.scalePhaseVolume('Aqueous', 0.1, 'm3') # corresponds to the initial porosity of 10%.
state_ic.scalePhaseVolume('Quartz', 0.882, 'm3')
state_ic.scalePhaseVolume('Calcite', 0.018, 'm3')

# Step 10: Scale the boundary condition state
state_bc.scaleVolume(1.0)

# Step 11: Create the mesh for the column
mesh = Mesh(ncells, xl, xr)

# Step 12: Create a chemical field object with every cell having state given by state_ic
field = ChemicalField(mesh.numCells(), state_ic)

# Step 13: Define the reactive transport modeling
rt = ReactiveTransportSolver(system)
rt.setMesh(mesh)
rt.setVelocity(v)
rt.setDiffusionCoeff(D)
rt.setBoundaryState(state_bc)
rt.setTimeStep(dt)
# Initialize reactive transport solver with chemical field object
rt.initialize(field)

# Step 14: Define the quantities that should be output for every cell, every time step
output = rt.output()
output.add('pH')
output.add('speciesMolality(H+)')
output.add('speciesMolality(Ca++)')
output.add('speciesMolality(Mg++)')
output.add('speciesMolality(HCO3-)')
output.add('speciesMolality(CO2(aq))')
output.add('phaseVolume(Calcite)')
output.add('phaseVolume(Dolomite)')

# Create the folder for the output text-files
os.system('mkdir -p results')
# Defined the name of the output files
output.filename('results/reactive-transport.txt')

# Step 15: Set initial time and counter of steps in time
t = 0.0
step = 0

# Reactive transport simulations in the cycle
while step <= nsteps:
    # Print the progress of the simulation
    print("Progress: {}/{} steps, {} min".format(step, nsteps, t/minute))
    
    # Perform one reactive transport time step
    rt.step(field)
    
    # Increment time step and number of time steps
    t += dt
    step += 1
