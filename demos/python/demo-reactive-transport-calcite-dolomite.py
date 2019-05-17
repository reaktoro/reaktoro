print('============================================================')
print('Make sure you have the following Python packages installed: ')
print('     numpy, matplotlib, joblib')
print('These can be installed with pip:')
print('     pip install numpy matplotlib joblib')
print('============================================================')

# Step 1: importing the required python packages
from reaktoro import *
from numpy import *
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
import os, time

# Step 2: Auxiliary time related constants
second = 1
minute = 60
hour = 60 * minute
day = 24 * hour
year = 365 * day

# Step 3: Parameters for the reactive transport simulation
xl = 0.0          # the x-coordinate of the left boundary
xr = 1.0          # the x-coordinate of the right boundary
ncells = 100      # the number of cells in the discretization
nsteps = 1200      # the number of steps in the reactive transport simulation
D  = 1.0e-9       # the diffusion coefficient (in units of m2/s)
v  = 1.0/day      # the fluid pore velocity (in units of m/s)
dx = (xr - xl)/ncells
dt = 5*minute     # the time step (in units of s)
T = 60.0 + 273.15  # the temperature (in units of K)
P = 100 * 1e5      # the pressure (in units of Pa)

alpha = v*dt/dx

dirichlet = False  # the parameter that determines whether Dirichlet BC must be used

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

# Step 7: Perform the reactive transport simulation
def simulate():

    # Step 7.1: Construct the chemical system with its phases and species
    db = Database('supcrt98.xml')

    editor = ChemicalEditor(db)

    #editor.addAqueousPhase(['H2O(l)', 'H+', 'OH-', 'Na+', 'Cl-', 'Ca++', 'Mg++', 'HCO3-', 'CO2(aq)', 'CO3--'])  # aqueous species are individually selected for performance reasons
    editor.addAqueousPhaseWithElementsOf('H2O NaCl CaCl2 MgCl2 CO2 SiO2 CaCO3')
    editor.addMineralPhase('Quartz')
    editor.addMineralPhase('Calcite')
    editor.addMineralPhase('Dolomite')

    # Step 7.2: Create the ChemicalSystem object using the configured editor
    system = ChemicalSystem(editor)

    # Step 7.3: Define the initial condition of the reactive transport modeling problem
    problem_ic = EquilibriumProblem(system)
    problem_ic.setTemperature(T)
    problem_ic.setPressure(P)
    problem_ic.add('H2O', 1.0, 'kg')
    problem_ic.add('NaCl', 0.7, 'mol')
    problem_ic.add('CaCO3', 10, 'mol')
    problem_ic.add('SiO2', 10, 'mol')

    # Step 7.4: Define the boundary condition of the reactive transport modeling problem
    problem_bc = EquilibriumProblem(system)
    problem_bc.setTemperature(T)
    problem_bc.setPressure(P)
    problem_bc.add('H2O', 1.0, 'kg')
    problem_bc.add('NaCl', 0.90, 'mol')
    problem_bc.add('MgCl2', 0.05, 'mol')
    problem_bc.add('CaCl2', 0.01, 'mol')
    problem_bc.add('CO2', 0.75, 'mol')

    options = EquilibriumOptions()
    # options.smart.reltol = 0.50
    # options.smart.abstol = 0.25
    options.smart.reltol = 0.1
    options.smart.abstol = 0.1

    # Step 7.5: Calculate the equilibrium states for the initial and boundary conditions
    state_ic = equilibrate(problem_ic)
    state_bc = equilibrate(problem_bc)

    # Step 7.6: Scale the volumes of the phases in the initial condition
    state_ic.scalePhaseVolume('Aqueous', 0.1, 'm3')
    state_ic.scalePhaseVolume('Quartz', 0.882, 'm3')
    state_ic.scalePhaseVolume('Calcite', 0.018, 'm3')

    # Scale the boundary condition state to 1 m3
    state_bc.scaleVolume(1.0, 'm3')


    # Step 7.7: Partitioning fluid and solid species

    # The number of elements in the chemical system
    nelems = system.numElements()

    # The indices of the fluid and solid species
    ifluid_species = system.indicesFluidSpecies()
    isolid_species = system.indicesSolidSpecies()

    # The concentrations of each element in each mesh cell (in the current time step)
    b = zeros((ncells, nelems))

    # The concentrations (mol/m3) of each element in the fluid partition, in each mesh cell
    bfluid = zeros((ncells, nelems))

    # The concentrations (mol/m3) of each element in the solid partition, in each mesh cell
    bsolid = zeros((ncells, nelems))

    # Initialize the concentrations (mol/m3) of the elements in each mesh cell
    b[:] = state_ic.elementAmounts()

    # Initialize the concentrations (mol/m3) of each element on the boundary
    b_bc = state_bc.elementAmounts()

    # Step 7.8: Create a list of chemical states for the mesh cells

    # The list of chemical states, one for each cell, initialized to state_ic
    states = [state_ic.clone() for _ in range(ncells)]

    # Step 7.9: Create the equilibrium solver object for the repeated equilibrium calculation

    # The interval [xl, xr] split into ncells
    x = linspace(xl, xr, ncells + 1)

    # The length of each mesh cell (in units of m)
    dx = (xr - xl)/ncells

    # Step 7.10: Create the equilibrium solver object for the repeated equilibrium calculation
    solver = EquilibriumSolver(system)

    # Step 7.11: The auxiliary function to create an output file each time step
    def outputstate():
        # Create the instance of ChemicalOutput class
        output = ChemicalOutput(system)

        # The number of digits in number of steps (e.g., 100 has 3 digits)
        ndigits = len(str(nsteps))

        # Provide the output file name, which will correspond
        output.filename('results/{}.txt'.format(str(step).zfill(ndigits)))

        # We define the columns' tags filled with the name of the quantities
        # The first column has a tag 'x' (which corresponds to the center coordinates of the cells )
        output.add('tag', 'x') # The value of the center coordinates of the cells

        # The rest of the columns correspond to the requested properties
        for quantity in output_quantities:
            output.add(quantity)

        # We update the file with states that correspond to the cells' coordinates stored in x
        output.open()
        for state, tag in zip(states, x):
            output.update(state, tag)
        output.close()

    # Step 7.12: Running the reactive transport simulation loop
    step = 0  # the current step number
    t = 0.0  # the current time (in seconds)

    while step <= nsteps:
        # Print the progress of the simulation
        print("Progress: {}/{} steps, {} min".format(step, nsteps, t/minute))

        # Collect the amounts of elements from fluid and solid partitions
        for icell in range(ncells):
            bfluid[icell] = states[icell].elementAmountsInSpecies(ifluid_species)
            bsolid[icell] = states[icell].elementAmountsInSpecies(isolid_species)

        # Transport each element in the fluid phase
        for j in range(nelems):
            #transport_fullimplicit(bfluid[:, j], dt, dx, v, D, b_bc[j])
            transport_implicit_explicit(bfluid[:, j], dt, dx, v, D, b_bc[j])

        # Update the amounts of elements in both fluid and solid partitions
        b[:] = bsolid + bfluid

        # Equilibrating all cells with the updated element amounts
        for icell in range(ncells):
            solver.solve(states[icell], T, P, b[icell])

        # Output the current state of the reactive transport calculation
        outputstate()

        # Increment time step and number of time steps
        t += dt
        step += 1

    print("Finished!")


# Step 10: Return a string for the title of a figure in the format Time: #h##m
def titlestr(t):
    t = t / minute   # Convert from seconds to minutes
    h = int(t) / 60  # The number of hours
    m = int(t) % 60  # The number of remaining minutes
    return 'Time: %2dh %2dm' % (h, m)

# Step 9: Generate figures for a result file
def plotfile(file):

    step = int(file.split('.')[0])

    print('Plotting figure', step, '...')

    t = step * dt
    filearray = loadtxt('results/' + file, skiprows=1)
    data = filearray.T

    ndigits = len(str(nsteps))

    plt.figure()
    plt.xlim(left=xl-0.02, right=xr+0.02)
    plt.ylim(bottom=2-0.1, top=9.0)
    plt.title(titlestr(t))
    plt.xlabel('Distance [m]')
    plt.ylabel('pH')
    plt.plot(data[0], data[1])
    plt.tight_layout()
    plt.savefig('figures/ph/{}.png'.format(str(step).zfill(ndigits)))

    plt.figure()
    plt.xlim(left=xl-0.02, right=xr+0.02)
    plt.ylim(bottom=-0.1, top=2.1)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.title(titlestr(t))
    plt.xlabel('Distance [m]')
    plt.ylabel('Mineral Volume [%$_{\mathsf{vol}}$]')
    plt.plot(data[0], data[7] * 100, label='Calcite')
    plt.plot(data[0], data[8] * 100, label='Dolomite')
    plt.legend(loc='center right')
    plt.tight_layout()
    plt.savefig('figures/calcite-dolomite/{}.png'.format(str(step).zfill(ndigits)))

    plt.figure()
    plt.yscale('log')
    plt.xlim(left=xl-0.02, right=xr+0.02)
    plt.ylim(bottom=0.5e-5, top=2)
    plt.title(titlestr(t))
    plt.xlabel('Distance [m]')
    plt.ylabel('Concentration [molal]')
    plt.plot(data[0], data[3], label='Ca++')
    plt.plot(data[0], data[4], label='Mg++')
    plt.plot(data[0], data[5], label='HCO3-')
    plt.plot(data[0], data[6], label='CO2(aq)')
    plt.plot(data[0], data[2], label='H+')
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig('figures/aqueous-species/{}.png'.format(str(step).zfill(ndigits)))

    plt.close('all')


# Step 8: Plot all result files and generate a video
def plot():
    # Plot all result files
    files = sorted(os.listdir('results'))
    #for file in files:
    #    plotfile(file);
    Parallel(n_jobs=16)(delayed(plotfile)(file) for file in files)

    # Create videos for the figures
    ffmpegstr = 'ffmpeg -y -r 30 -i figures/{0}/%03d.png -codec:v mpeg4 -flags:v +qscale -global_quality:v 0 videos/{0}.mp4'
    os.system(ffmpegstr.format('calcite-dolomite'))
    os.system(ffmpegstr.format('aqueous-species'))
    os.system(ffmpegstr.format('ph'))

# Step 7.10.2: Solve a tridiagonal matrix equation using Thomas algorithm (or TriDiagonal Matrix Algorithm (TDMA))
def thomas(a, b, c, d):
    n = len(d)
    c[0] /= b[0]
    for i in range(1, n - 1):
        c[i] /= b[i] - a[i]*c[i - 1]
    d[0] /= b[0]
    for i in range(1, n):
        d[i] = (d[i] - a[i]*d[i - 1])/(b[i] - a[i]*c[i - 1])
    x = d
    for i in reversed(range(0, n - 1)):
        x[i] -= c[i]*x[i + 1]
    return x

# Step 7.10.1: Perform a transport step
# Both diffusion and convection are treated implicitly
def transport_fullimplicit(u, dt, dx, v, D, ul):

    # Number of DOFs
    n = len(u)
    alpha = D*dt/dx**2
    beta = v*dt/dx

    # Upwind finite-volume scheme
    a = full(n, -beta - alpha)
    b = full(n, 1 + beta + 2*alpha)
    c = full(n, -alpha)

    # Set the boundary condition on the left cell
    if dirichlet:
        # Use Dirichlet BC boundary conditions
        b[0] = 1.0
        c[0] = 0.0
        u[0] = ul

    else:
        # Flux boundary conditions (implicit scheme for the advection)
        # Left boundary
        b[0] = 1 + alpha + beta
        c[0] = -alpha # stays the same as it is defined -alpha
        u[0] += beta * ul # = dt/dx * v * g, flux that we prescribe is equal v * ul

    # Right boundary is free
    a[-1] = - beta
    b[-1] = 1 + beta

    # Solve a tridiagonal matrix equation
    thomas(a, b, c, u)

# Step 7.10.1: Perform a transport step
# Both diffusion and convection are treated implicitly
def transport_fullimplicit_left_zeroflux(u, dt, dx, v, D, ul):

    # Number of DOFs
    n = len(u)
    alpha = D*dt/dx**2
    beta = v*dt/dx

    # Upwind finite-volume scheme
    a = full(n, -beta - alpha)
    b = full(n, 1 + beta + 2*alpha)
    c = full(n, -alpha)

    # Set the boundary condition on the left cell
    if dirichlet:
        # Use Dirichlet BC boundary conditions
        b[0] = 1.0
        c[0] = 0.0
        u[0] = ul
    else:
        # Flux boundary conditions (implicit scheme for the advection)
        # Left boundary
        u[0] += beta * ul
        b[0] = 1 + alpha + beta
        c[0] = -alpha # stays the same as it is defined -alpha

    # Right boundary
    a[-1] = - beta - alpha
    b[-1] = 1 + alpha

    # Solve a tridiagonal matrix equation
    thomas(a, b, c, u)

# Diffusion is treated implicitly and convection explicitly
def transport_implicit_explicit(u, dt, dx, v, D, ul):

    # Number of DOFs
    n = len(u)
    alpha = D*dt/dx**2
    beta = v*dt/dx

    # Upwind finite-volume scheme
    a = full(n, - alpha)
    b = full(n, 1 + 2*alpha)
    c = full(n, - alpha)

    u0 = u

    # Set the boundary condition on the left cell
    if dirichlet:
        # Use Dirichlet BC boundary conditions
        b[0] = 1.0
        c[0] = 0.0
        u[0] = ul

    else:
        # Use flux boundary conditions (explicit scheme for the advection)
        # Left boundary
        b[0] = 1 + alpha
        c[0] = -alpha
        u[0] = (1 - beta) * u0[0] + beta * ul

    # Contribution of the explicit advection to the RHS
    for i in range(1, n-1):
        u[i] = beta * u0[i-1] + (1 + beta) * u0[i]

    # Right boundary
    a[-1] = 0
    b[-1] = 1
    u[-1] = beta*u0[-2] + (1 - beta)*u0[-1]

    # Solve a tridiagonal matrix equation
    thomas(a, b, c, u)

# Step 6: Creating folders for the results' output
def make_results_folders():
    os.system('mkdir -p results')
    os.system('mkdir -p figures/ph')
    os.system('mkdir -p figures/aqueous-species')
    os.system('mkdir -p figures/calcite-dolomite')
    os.system('mkdir -p videos')

# Step 5: Define the main function
if __name__ == '__main__':

    # Create folders for the results
    make_results_folders()

    # Run the reactive transport simulations
    simulate()

    # Plotting the result
    plot()
