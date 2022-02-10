print('============================================================')
print('Make sure you have the following Python packages installed: ')
print('     numpy, matplotlib, joblib')
print('These can be installed with pip:')
print('     pip install numpy matplotlib joblib')
print('============================================================')

# Step 1: importing the required python packages
from reaktoro import *
from numpy import *
import plotting
import os, time

# Step 2: Auxiliary time related constants
second = 1
minute = 60
hour = 60 * minute
day = 24 * hour
week = 7 * day
year = 365 * day

# Step 3: Parameters for the reactive transport simulation
xl = 0.0               # the x-coordinate of the left boundary
xr = 1.0               # the x-coordinate of the right boundary
ncells = 100           # the number of cells in the discretization
nsteps = 300            # the number of steps in the reactive transport simulation
D  = 1.0e-9            # the diffusion coefficient (in units of m2/s)
v  = 1.0/week          # the fluid pore velocity (in units of m/s)
dx = (xr - xl)/ncells  # the length of the mesh cells
dt = 30*minute         # the time step (30 minutes in in units of s)
T = 25.0               # the temperature (in units of degC)
P = 100                # the pressure (in units of bar)
rsa = 1.0

alpha = v*dt/dx

dirichlet = False  # the parameter that determines whether Dirichlet BC must be used

# set which kinetic species are active
include_dolomite = True
include_calcite = True

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
output_quantities_map = {"pH" : 1,
                         "speciesMolality(H+)" : 2,
                         "speciesMolality(Ca++)" : 3,
                         "speciesMolality(Mg++)" : 4,
                         "speciesMolality(HCO3-)" : 5,
                         "speciesMolality(CO2(aq))" : 6,
                         "phaseVolume(Calcite)" : 7,
                         "phaseVolume(Dolomite)" : 8}

folder = 'results-rt-kinetic-dolomitization'

# Step 7: Perform the reactive transport simulation
def simulate():

    # Step 7.1: Construct the chemical system with its phases and species
    db = Database('supcrt98.xml')

    editor = ChemicalEditor(db)

    # Define phases and corresponding to these phases species
    editor.addAqueousPhaseWithElementsOf("H2O NaCl CaCO3 MgCO3")
    editor.addGaseousPhase(["H2O(g)", "CO2(g)"])
    editor.addMineralPhase("Calcite")
    editor.addMineralPhase("Dolomite")
    editor.addMineralPhase("Halite")

    # define chemical system
    system = ChemicalSystem(editor)

    # Partition system to kinetic species
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

    # Create the ChemicalSystem and ReactionSystem objects using the configured editor
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

    # Step 7.3: Define the initial condition of the reactive transport modeling problem
    problem_ic = EquilibriumProblem(system)
    problem_ic.setTemperature(T, 'celsius')
    problem_ic.setPressure(P, 'bar')
    problem_ic.add('H2O', 1.0, 'kg')
    problem_ic.add('NaCl', 0.1, 'mol')
    problem_ic.add('CO2', 0.001, 'mol')
    problem_ic.add('Calcite', 10, 'kg')
    problem_ic.add('Dolomite', 10, 'kg')

    # Step 7.4: Define the boundary condition of the reactive transport modeling problem
    problem_bc = EquilibriumProblem(system)
    problem_bc.setTemperature(T, 'celsius')
    problem_bc.setPressure(P, 'bar')
    problem_bc.add('H2O', 1.0, 'kg')
    problem_bc.add('NaCl', 0.90, 'mol')
    problem_bc.add('MgCl2', 0.05, 'mol')
    problem_bc.add('CaCl2', 0.01, 'mol')
    problem_bc.add('CO2', 0.75, 'mol')

    # Step 7.5: Calculate the equilibrium states for the initial and boundary conditions
    state_ic = equilibrate(problem_ic)
    state_bc = equilibrate(problem_bc)

    # Step 7.6: Scale the volumes of the phases in the initial condition
    state_ic.scalePhaseVolume('Aqueous', 0.1, 'm3')
    state_ic.scalePhaseVolume('Calcite', 0.45, 'm3')
    state_ic.scalePhaseVolume('Dolomite', 0.45, 'm3')

    # Scale the boundary condition state to 1 m3
    state_bc.scaleVolume(1.0, 'm3')

    # Step 7.7: Partitioning fluid and solid species

    # The number of elements in the chemical system
    nelems = system.numElements()

    # The indices of the fluid and solid species
    ifs = partition.indicesFluidSpecies()
    iss = partition.indicesSolidSpecies()

    # The concentrations of each element in each mesh cell (in the current time step)
    be = zeros((ncells, nelems))

    # The concentrations (mol/m3) of each element in the fluid partition, in each mesh cell
    bef = zeros((ncells, nelems))

    # The concentrations (mol/m3) of each element in the solid partition, in each mesh cell
    bes = zeros((ncells, nelems))

    # # Initialize the concentrations (mol/m3) of the elements in each mesh cell
    be[:] = state_ic.elementAmounts()

    # Initialize the concentrations (mol/m3) of each element on the boundary
    be_bc = state_bc.elementAmountsInSpecies(ifs)

    # Step 7.8: Create a list of chemical states for the mesh cells
    # The list of chemical states, one for each cell, initialized to state_ic
    states = [state_ic.clone() for _ in range(ncells)]

    # Step 7.9: Create the equilibrium solver object for the repeated equilibrium calculation

    # The interval [xl, xr] split into ncells
    x = linspace(xl, xr, ncells + 1)

    # The length of each mesh cell (in units of m)
    dx = (xr - xl)/ncells

    # Step 7.10: Create the equilibrium solver object for the repeated equilibrium calculation
    solver = KineticSolver(reactions, partition)

    # Step 7.11: The auxiliary function to create an output file each time step
    def outputstate():
        # Create the instance of ChemicalOutput class
        output = ChemicalOutput(system)

        # The number of digits in number of steps (e.g., 100 has 3 digits)
        ndigits = len(str(nsteps))

        # Provide the output file name, which will correspond
        output.filename(folder + '/step-{}.txt'.format(str(step).zfill(ndigits)))

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
    outputstate()

    while step <= nsteps:
        # Print the progress of the simulation
        print("Progress: {}/{} steps, {} min".format(step, nsteps, t/minute))

        # Collect the amounts of elements from fluid and solid partitions
        for icell in range(ncells):

            bef[icell] = states[icell].elementAmountsInSpecies(ifs)
            bes[icell] = states[icell].elementAmountsInSpecies(iss)

        # Get the porosity of the left boundary cell
        bc_cell = 0
        phi_bc = states[bc_cell].properties().fluidVolume().val / states[bc_cell].properties().volume().val

        # Transport each element in the fluid phase
        for j in range(nelems):
            transport_fullimplicit(bef[:, j], dt, dx, v, D, phi_bc * be_bc[j])

        # Update the amounts of elements in both fluid and solid partitions
        be[:] = bes + bef

        # Equilibrating all cells with the updated element amounts
        for icell in range(ncells):
            solver.solve(states[icell], t, dt, be[icell])

        # Increment time step and number of time steps
        t += dt
        step += 1

        # Output the current state of the reactive transport calculation
        outputstate()

    print("Finished!")

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
    os.system('mkdir -p ' + folder)
    os.system('mkdir -p figures-' + folder + '/ph')
    os.system('mkdir -p figures-' + folder + '/aqueous-species')
    os.system('mkdir -p figures-' + folder + '/calcite-dolomite')
    os.system('mkdir -p videos-' + folder)

x = linspace(xl, xr, ncells)

# Step 5: Define the main function
if __name__ == '__main__':

    # Create folders for the results
    make_results_folders()

    # Run the reactive transport simulations
    simulate()

    # Plotting the result
    plotting.plot_rt_kinetic_dolomitization(folder, x, nsteps, dt, output_quantities_map)
