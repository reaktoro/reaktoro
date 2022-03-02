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

# Step 1: Import the reaktoro Python package (and other packages)
from reaktoro import *
from numpy import *
import os
import plotting

#------------------------------------------------------------------------------#
# Problems parameters
#------------------------------------------------------------------------------#

# Step 2: Initialise auxiliary time-related constants
second = 1
minute = 60
hour = 60 * minute
day = 24 * hour
week = 7 * day
year = 365 * day

# Step 3: Define parameters for the reactive transport simulation
xl = 0.0          # the x-coordinate of the left boundary
xr = 1.0          # the x-coordinate of the right boundary
ncells = 100      # the number of cells in the discretization
nsteps = 300      # the number of steps in the reactive transport simulation
D  = 1.0e-9       # the diffusion coefficient (in units of m2/s)
v  = 1.0/week     # the fluid pore velocity (1 m/week in units of m/s)
dt = 30*minute    # the time step (30 minutes in units of s)
T = 25.0          # the temperature (in units of degC)
P = 100           # the pressure (in units of bar)
rsa = 1.0         # specific surface area

# set which kinetic species are active
include_dolomite = True
include_calcite = True

dx = (xr - xl)/ncells
alpha = v*dt/dx
ndigits = len(str(nsteps))
folder = 'results-rtsolver-kinetic-dolomitization-with-temperature-change'

# Step 4: The list of quantities to be output for each mesh cell, each time step
output_quantities = """
    pH
    speciesMolality(H+)
    speciesMolality(Ca++)
    speciesMolality(Mg++)
    speciesMolality(HCO3-)
    speciesMolality(CO2(aq))
    phaseVolume(Calcite)
    phaseVolume(Dolomite)
""".split()
output_quantities_map = {"pH" : 0,
                         "speciesMolality(H+)" : 1,
                         "speciesMolality(Ca++)" : 2,
                         "speciesMolality(Mg++)" : 3,
                         "speciesMolality(HCO3-)" : 4,
                         "speciesMolality(CO2(aq))" : 5,
                         "phaseVolume(Calcite)" : 6,
                         "phaseVolume(Dolomite)" : 7}
#------------------------------------------------------------------------------#
# Auxiliary functions
#------------------------------------------------------------------------------#

# Step 6: Creating folders for the results' output
def make_results_folders():
    os.system('mkdir -p ' + folder)
    os.system('mkdir -p figures-' + folder + '/ph')
    os.system('mkdir -p figures-' + folder + '/aqueous-species')
    os.system('mkdir -p figures-' + folder + '/calcite-dolomite')
    os.system('mkdir -p videos-' + folder)

def construct_temperature_sequence(T_bc, T_init, temperature_steps):

    # construct array with log temperature distribution
    T0 = np.logspace(T_bc, T_init, ncells + 1, endpoint = True, base =  1.7)
    T0 = T_init - np.flip(T0/np.max(T0) * (T_init-T_bc))
    T0 = np.round(T0,2)

    # replicate for every timestep
    T = np.tile(T0, (nsteps+1, 1))

    # indices at which to change temperature
    inds = np.round((ncells+1)/temperature_steps,0)

    # propagate temperature
    for i in range(1, temperature_steps):
        T[int(inds*i):, 1:] = T[int(inds*i):, :-1]

    # restore initial temperature in first step
    T[0,:] = T_init

    return T

#------------------------------------------------------------------------------#
# Reactive transport simulations
#------------------------------------------------------------------------------#

temperature_steps = 20
T_init = 140 + 273.15  # temperature (in units of K)
T_bc = 50 + 273.15

# Define temperature array
temperatures = construct_temperature_sequence(T_bc, T_init, temperature_steps)

# Step 4: Construct the chemical system with its phases and species
db = Database('supcrt98.xml')
editor = ChemicalEditor(db)
editor.addAqueousPhaseWithElementsOf("H2O NaCl CaCO3 MgCO3")
editor.addGaseousPhase(["H2O(g)", "CO2(g)"])
editor.addMineralPhase("Calcite")
editor.addMineralPhase("Dolomite")
editor.addMineralPhase("Halite")

editor.addMineralReaction("Calcite") \
    .setEquation("Calcite = Ca++ + CO3--") \
    .addMechanism("logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol") \
    .addMechanism("logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0") \
    .setSpecificSurfaceArea(rsa, "cm2/g")

editor.addMineralReaction("Dolomite") \
    .setEquation("Dolomite = Ca++ + Mg++ + 2*CO3--") \
    .addMechanism("logk = -7.53 mol/(m2*s); Ea = 52.2 kJ/mol") \
    .addMechanism("logk = -3.19 mol/(m2*s); Ea = 36.1 kJ/mol; a[H+] = 0.5") \
    .setSpecificSurfaceArea(rsa, "cm2/g")

# Step 5: Create the ChemicalSystem and ReactionSystem objects using the configured editor
system = ChemicalSystem(editor)
reactions = ReactionSystem(editor)

partition = Partition(system)
if (include_dolomite == True) & (include_calcite == False):
    partition.setKineticSpecies(["Dolomite"])
    print("only Dolomite")
elif (include_dolomite == False) & (include_calcite == True):
    partition.setKineticSpecies(["Calcite"])
    print("only Calcite")
elif (include_dolomite == True) & (include_calcite == True):
    partition.setKineticSpecies(["Calcite", "Dolomite"])
    print("Calcite & Dolomite")

# Step 6: Define the initial condition of the reactive transport modeling problem
problem_ic = EquilibriumProblem(system)
problem_ic.setTemperature(T_init, 'kelvin')
problem_ic.setPressure(P, 'bar')
problem_ic.add('H2O', 1.0, 'kg')
problem_ic.add('NaCl', 0.1, 'mol')
problem_ic.add('CO2', 0.001, 'mol')
problem_ic.add('Calcite', 10, 'kg')
problem_ic.add('Dolomite', 10, 'kg')

# Step 7: Define the boundary condition of the reactive transport modeling problem
problem_bc = EquilibriumProblem(system)
problem_bc.setTemperature(T_bc, 'kelvin')
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
state_ic.scalePhaseVolume('Aqueous', 0.1, 'm3')
state_ic.scalePhaseVolume('Calcite', 0.45, 'm3')
state_ic.scalePhaseVolume('Dolomite', 0.45, 'm3')

# Step 10: Scale the boundary condition state
state_bc.scaleVolume(1.0, 'm3')

# Step 11: Create the mesh for the column
mesh = Mesh(ncells + 1, xl, xr)

# The interval [xl, xr] split into ncells
x = linspace(xl, xr, ncells + 1)

# Step 12: Create a chemical field object with every cell having state given by state_ic
field = ChemicalField(mesh.numCells(), state_ic)

# Step 13: Initialize the reactive transport solver
rt = ReactiveTransportSolver(reactions, partition)
rt.setMesh(mesh)
rt.setVelocity(v)
rt.setDiffusionCoeff(D)
rt.setBoundaryState(state_bc)
rt.setTimeStep(dt)
rt.initialize()

# Step 14: Set the output of the reactive transport simulation
output = rt.output()
output.add("pH")
output.add("speciesMolality(H+)")
output.add("speciesMolality(Ca++)")
output.add("speciesMolality(Mg++)")
output.add("speciesMolality(HCO3-)")
output.add("speciesMolality(CO2(aq))")
output.add("phaseVolume(Calcite)")
output.add("phaseVolume(Dolomite)")
output.filename(folder + '/step.txt')  # Set the name of the output files

make_results_folders()

# Step 15: Perform given number of reactive tranport steps
t = 0.0  # current time variable
step = 0  # current number of steps

while step <= nsteps:  # step until the number of steps are achieved
    # Print the progress of the simulation
    print("Progress: {}/{} steps, {} min".format(step, nsteps, t/minute))

    # Update the array of temperatures in the chemical field on the next step
    if step > 0:
        for icell in range(ncells + 1):
            field[icell].setTemperature(temperatures[step, icell])

    # Perform one reactive transport time step
    rt.step(field)

    # Increment time step and number of time steps
    t += dt
    step += 1

plotting.plot_rt_kinetic_dolomitization(folder, x, nsteps, dt, output_quantities_map)
