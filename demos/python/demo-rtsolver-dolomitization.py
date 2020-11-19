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
from joblib import Parallel, delayed
import matplotlib.pyplot as plt

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
nsteps = 600      # the number of steps in the reactive transport simulation
D  = 1.0e-9       # the diffusion coefficient (in units of m2/s)
v  = 1.0/week     # the fluid pore velocity (1 m/week in units of m/s)
dt = 30*minute    # the time step (30 minutes in units of s)
T = 60.0          # the temperature (in units of degC)
P = 100           # the pressure (in units of bar)

dx = (xr - xl)/ncells
alpha = v*dt/dx
ndigits = len(str(nsteps))
folder = 'results-rtsolver-dolomitization'

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

#------------------------------------------------------------------------------#
# Auxiliary functions
#------------------------------------------------------------------------------#

# Step 10: Return a string for the title of a figure in the format Time: #h##m
def titlestr(t):
    t = t / minute   # Convert from seconds to minutes
    h = int(t) / 60  # The number of hours
    m = int(t) % 60  # The number of remaining minutes
    return 'Time: {:>3}h{:>2}m'.format(h, str(m).zfill(2))

# Step 6: Creating folders for the results' output
def make_results_folders():
    os.system('mkdir -p results')
    os.system('mkdir -p figures/ph')
    os.system('mkdir -p figures/aqueous-species')
    os.system('mkdir -p figures/calcite-dolomite')
    os.system('mkdir -p videos')

# Step 8: Plot all result files and generate a video
def plot():
    # Plot all result files
    files = sorted(os.listdir('results'))
    Parallel(n_jobs=16)(delayed(plotfile)(file) for file in files)
    # Create videos for the figures
    ffmpegstr = 'ffmpeg -y -r 30 -i figures/{0}/%0' + str(ndigits) + 'd.png -codec:v mpeg4 -flags:v +qscale -global_quality:v 0 videos/{0}.mp4'
    os.system(ffmpegstr.format('calcite-dolomite'))
    os.system(ffmpegstr.format('aqueous-species'))
    os.system(ffmpegstr.format('ph'))

def plotfile(file):

    step = int((file.split('.')[0]).split('-')[1])

    print('Plotting figure', step, '...')

    t = step * dt
    filearray = loadtxt('results/' + file, skiprows=1)
    data = filearray.T

    ndigits = len(str(nsteps))

    x = linspace(xl, xr, ncells)

    plt.figure()
    plt.xlim(left=-0.02, right=0.52)
    plt.ylim(bottom=2.5, top=10.5)
    plt.title(titlestr(t))
    plt.xlabel('Distance [m]')
    plt.ylabel('pH')
    plt.plot(x, data[0])
    plt.tight_layout()
    plt.savefig('figures/ph/{}.png'.format(str(step).zfill(ndigits)))

    plt.figure()
    plt.xlim(left=-0.02, right=0.52)
    plt.ylim(bottom=-0.1, top=2.1)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.title(titlestr(t))
    plt.xlabel('Distance [m]')
    plt.ylabel('Mineral Volume [%$_{\mathsf{vol}}$]')
    plt.plot(x, data[6] * 100, label='Calcite')
    plt.plot(x, data[7] * 100, label='Dolomite')
    plt.legend(loc='center right')
    plt.tight_layout()
    plt.savefig('figures/calcite-dolomite/{}.png'.format(str(step).zfill(ndigits)))

    plt.figure()
    plt.yscale('log')
    plt.xlim(left=-0.02, right=0.52)
    plt.ylim(bottom=0.5e-5, top=2)
    plt.title(titlestr(t))
    plt.xlabel('Distance [m]')
    plt.ylabel('Concentration [molal]')
    plt.plot(x, data[2], label='Ca++')
    plt.plot(x, data[3], label='Mg++')
    plt.plot(x, data[4], label='HCO3-')
    plt.plot(x, data[5], label='CO2(aq)')
    plt.plot(x, data[1], label='H+')
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig('figures/aqueous-species/{}.png'.format(str(step).zfill(ndigits)))

    plt.close('all')

#------------------------------------------------------------------------------#
# Reactive transport simulations
#------------------------------------------------------------------------------#

# Step 4: Construct the chemical system with its phases and species
db = Database('supcrt98.xml')
editor = ChemicalEditor(db)
editor.addAqueousPhase('H2O(l) H+ OH- Na+ Cl- Ca++ Mg++ HCO3- CO2(aq) CO3--')
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
state_bc.scaleVolume(1.0, 'm3')

# Step 11: Create the mesh for the column
mesh = Mesh(ncells, xl, xr)

# Step 12: Create a chemical field object with every cell having state given by state_ic
field = ChemicalField(mesh.numCells(), state_ic)

# Step 13: Initialize the reactive transport solver
rt = ReactiveTransportSolver(system)
rt.setMesh(mesh)
rt.setVelocity(v)
rt.setDiffusionCoeff(D)
rt.setBoundaryState(state_bc)
rt.setTimeStep(dt)
rt.initialize(field)

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
output.filename('results/rtsolver.txt')  # Set the name of the output files

make_results_folders()

# Step 15: Perform given number of reactive tranport steps
t = 0.0  # current time variable
step = 0  # current number of steps

while step <= nsteps:  # step until the number of steps are achieved
    # Print the progress of the simulation
    print("Progress: {}/{} steps, {} min".format(step, nsteps, t/minute))

    # Perform one reactive transport time step
    rt.step(field)

    # Increment time step and number of time steps
    t += dt
    step += 1


plot()
