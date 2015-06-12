import math
import numpy
from dolfin import *
from reaktoro.core import *
from reaktoro.core.mobility import Mobility
from reaktoro.transport.transport import TransportSolver
from reaktoro.core.porosity import Porosity
from numpy import empty
import time

def elementAmounts(states, result):
    for k in range(len(states)):
        result[k] = states[k].elementAmounts()

def elementAmountsInSpecies(states, ispecies, result):
    for k in range(len(states)):
        bspecies = states[k].elementAmountsInSpecies(ispecies)
        result[k] = bspecies.array()

def phaseVolumes(state):
    T = state.temperature()
    P = state.presure()
    n = state.speciesAmounts()
    return state.system().phaseVolumes(T, P, n)

def porosity(state, mobility):
    v = phaseVolumes(state)
    vs = sum([v[i] for i in mobility.indicesSolidPhases()])
    return 1.0 - vs

class ChemicalDirichletBC:
    def __init__(self, function_space, state, mobility, boundary):
        self.fluid = mobility.fluid(state)
        self.values = Function(function_space)
        self.dirichlet = DirichletBC(function_space, 1, boundary)
        self.dirichlet.apply(self.values.vector())
        self.dofs = numpy.where(self.values.vector() == 1)[0]
        self.dirichlet = DirichletBC(function_space, self.values, boundary)

    def elementDirichletBC(self, ielement, ispecies, porosity):
        bval = self.fluid.elementAmountInSpecies(ielement, ispecies)
        phi = porosity.vector()[self.dofs]
        self.values.vector()[self.dofs] = phi*bval
        return self.dirichlet

    def speciesDirichletBC(self, ispecies, porosity):
        nval = self.fluid.speciesAmount(ispecies)
        phi = porosity.vector()[self.dofs]
        self.values.vector()[self.dofs] = phi*nval
        return self.dirichlet


class ChemicalTransportResult(object):
    class EquilibriumResult(object):
        def __init__(self):
            # The number of iterations for equilibrium calculations performed at every degree of freedom
            self.iterations = []

            # The elapsed seconds for equilibrium calculations performed at every degree of freedom
            self.seconds = []

    class KineticsResult(object):
        def __init__(self):
            # The number of time-steps for kinetic calculations performed at every degree of freedom
            self.timesteps = []

            # The elapsed seconds for kinetic calculations performed at every degree of freedom
            self.seconds = []

    def __init__(self):
        # The flag that indicates if the chemical transport calculations succeeded
        self.succeeded = False

        # The total time in seconds of the chemical transport calculation
        self.seconds = 0

        # The total time in seconds spent on equilibrium calculations
        self.seconds_equilibrium = 0

        # The total time in seconds spent on finite element assembly operations
        self.seconds_assemble = 0

        # The total time in seconds spent on linear systems
        self.seconds_linear_systems = 0

        # The result of the equilibrium calculations
        self.equilibrium = ChemicalTransportResult.EquilibriumResult()

        # The result of the kinetic calculations
        self.kinetics = ChemicalTransportResult.KineticsResult()


class _ChemicalTransportSolver(object):
    def __init__(self, system, mobility):
        self.system = system
        self.mobility = mobility
        self.partition = Partition(system)
        self.num_species  = system.numSpecies()
        self.num_elements = system.numElements()
        self.num_fluid_phases = len(mobility.indicesFluidPhases())
        self.num_solid_phases = len(mobility.indicesSolidPhases())
        self.velocity = [Constant(0.0) for i in xrange(self.num_fluid_phases)]
        self.diffusion = [Constant(0.0) for i in xrange(self.num_fluid_phases)]
        self.source = [Constant(0.0) for i in xrange(self.num_fluid_phases)]
        self.boundary_conditions = []
        self.initialized = False

    def setVelocity(self, velocity):
        if self.num_fluid_phases > 1:
            if type(velocity) is not list:
                raise TypeError('Expecting a list of velocity fields for each fluid phase.')
            assert len(velocity) == self.num_fluid_phases, \
                'There are %d fluid phases, but only %d velocity fields given.' % \
                    (self.num_fluid_phases, len(velocity))
        else:
            self.velocity = velocity
        self.initialized = False

    def setDiffusion(self, diffusion):
        self.diffusion = diffusion
        self.initialized = False

    def setSource(self, source):
        self.source = source
        self.initialized = False

    def setPartition(self, partition):
        self.partition = partition

    def addBoundaryCondition(self, state, boundary):
        self.boundary_conditions.append((state, boundary))

    def setTemperatures(self, temperatures):
        for k in range(self.num_dofs):
            self.states[k].setTemperature(temperatures[k])

    def setPressures(self, pressures):
        for k in range(self.num_dofs):
            self.states[k].setPressure(pressures[k])

    def initialize(self, field):
        # Initialise the function space used in the ChemicalField instance
        self.V = field.functionSpace()

        # Initialize the list of global `dof` indices in the current process
        self.dofs = self.V.dofmap().dofs()

        # The number of degrees-of-freedom in the current process
        self.num_dofs = len(self.dofs)

        # Check if the user provided any boundary conditions
        if self.boundary_conditions is []:
            RuntimeError('Failed to initialize ChemicalTransportSolver. \
            No boundary conditions have been provided.')

        # Create the ChemicalDirichletBC instances
        self.bcs = [ChemicalDirichletBC(self.V, state, self.mobility, boundary) \
                    for (state, boundary) in self.boundary_conditions]

        # The auxiliary Function instance used for transport steps
        self.u = Function(self.V)

        # The auxiliary Function instance used for outputting
        self.output = Function(self.V)

        # The chemical equilibrium solver
        self.equilibrium = EquilibriumSolver(system)
        self.equilibrium.setPartition(self.partition)

        # Initialise the indices of the equilibrium and kinetic species
        self.ispecies_e  = self.partition.indicesEquilibriumSpecies()
        self.ispecies_k  = self.partition.indicesKineticSpecies()

        # Initialise the indices of the fluid and solid species
        self.ispecies_f  = self.mobility.indicesFluidSpeciesInEachFluidPhase()
        self.ispecies_s  = self.mobility.indicesSolidSpecies()

        # Initialise the indices of the equilibrium-fluid and equilibrium-solid species
        self.ispecies_ef = [sorted(set(self.ispecies_e) & set(indices)) for indices in self.ispecies_f]
        self.ispecies_es = sorted(set(self.ispecies_e) & set(self.ispecies_s))

        # Initialise the indices of the kinetic-fluid and kinetic-solid species
        self.ispecies_kf = [sorted(set(self.ispecies_k) & set(indices)) for indices in self.ispecies_f]
        self.ispecies_ks = sorted(set(self.ispecies_k) & set(self.ispecies_s))

        # Initialise the indices of fluid and solid phases
        self.iphases_f = mobility.indicesFluidPhases()
        self.iphases_s = mobility.indicesSolidPhases()

        # Initialise the number of fluid and solid phases
        self.num_fluid_phases = len(self.iphases_f)
        self.num_solid_phases = len(self.iphases_s)

        # Initialise the arrays of molar amounts of each element for each degree of freedom
        self.be  = numpy.zeros((self.num_elements, self.num_dofs)) # in the equilibrium partition
        self.bef = numpy.zeros((self.num_elements, self.num_dofs)) # in the equilibrium-fluid partition
        self.bes = numpy.zeros((self.num_elements, self.num_dofs)) # in the equilibrium-solid partition

        # Initialise the dolfin Function instances for the saturation field of each fluid phase
        self.saturation = [Function(self.V) for i in self.iphases_f]

        # Initialise the numpy array for the saturation field of each fluid phase
        self.saturation_array = numpy.zeros((self.num_fluid_phases, self.num_dofs))

        # Initialise the dolfin Function instance for the porosity field
        self.porosity = Function(self.V)

        # Initialise the numpy array for the porosity field
        self.porosity_array = numpy.zeros(self.num_dofs)

        # Define the pore velocity of each fluid phase
        pore_velocities = [self.velocity/(self.porosity*s) for s in self.saturation]

        # Initialise the transport solver
        self.transport = TransportSolver()
        self.transport.setVelocity(pore_velocity)
        self.transport.setDiffusion(self.diffusion)

        # Initialise the ChemicalTransportResult instance
        self.result = ChemicalTransportResult()
        self.result.equilibrium.iterations = Function(self.V)
        self.result.equilibrium.seconds = Function(self.V)
        self.result.kinetics.timesteps = Function(self.V)
        self.result.kinetics.seconds = Function(self.V)
        self.result.equilibrium._iterations = numpy.empty(self.num_dofs)
        self.result.equilibrium._seconds = numpy.empty(self.num_dofs)

    def updatePorositySaturation(self, field):
        # Calculate the porosity and saturation of each fluid phase in every degree-of-freedom
        for i in range(self.num_dofs):
            state = self.states[i]
            T = state.temperature()
            P = state.pressure()
            n = state.speciesAmounts()
            v = system.phaseVolumes(T, P, n).val()
            volume_fluid = sum([v[i] for i in self.iphases_f])
            volume_solid = sum([v[i] for i in self.iphases_s])
            self.porosity_array[k] = 1.0 - volume_solid
            self.saturation_array[:, k] = [v[i]/volume_fluid for i in self.iphases_f]

        # Update the Function instance corresponding to the porosity field
        self.porosity.vector()[:] = self.porosity_array

        # Update the Function instances corresponding to the saturation fields
        for i in range(self.num_fluid_phases):
            self.saturation[i].vector()[:] = self.saturation_array[i]


    def transportEquilibriumElementInFluidPhase(self, ielement, iphase, field, dt):
        # Create the Dirichlet boundary condition for the current element in the equilibrium-fluid partition
        bc_element = self.bc.elementDirichletBC(ielement, self.ispecies_ef, self.phi)

        # Set the boundary condition in the tranport solver
        self.transport.setBoundaryCondition(bc_element)

        # Set initial condition for the transport equation of the current element
        self.u.vector()[:] = self.bef[ielement]

        # Transport the current element of the equilibrium-fluid partition
        self.transport.step(self.u, dt, bc_element)

        # Extract the result to the array of element amounts `bef`
        self.bef[ielement] = self.u.vector().get_local()

    def transportEquilibriumFluidSpecies(self, field, dt):
        # Calculate the amounts of elements in the equilibrium-fluid partition
        field.elementAmountsInSpecies(self.ispecies_ef, self.bef)

        # Calculate the amounts of elements in the equilibrium-solid partition
        field.elementAmountsInSpecies(self.ispecies_es, self.bes)

        # Iterate over all elements in the equilibrium-fluid partition and tranport them
        for j in range(self.num_elements):
            self.transportEquilibriumElementInFluidPhase(j, field, dt)

        # Calculate the new amounts of elements in the equilibrium partition
        numpy.add(self.bes, self.bef, self.be)

        # Equilibrate the chemical states for every degree-of-freedom
        self.equilibrate(field)

    def equilibrate(self, field):
        tbegin = time.time()
        states = field.states()
        iterations = self.result.equilibrium._iterations
        seconds = self.result.equilibrium._seconds
        for i in range(len(states)):
            # Perform the equilibrium calculation at current degree of freedom
            result = self.equilibrium.solve(states[i], self.be[i])

            # Check if the equilibrium calculation was successful
            if not result.optimum.succeeded:
                raise RuntimeError("Failed to calculate equilibrium state at dof (%d), under \
                temperature %f K, pressure %f Pa, and element molar amounts %s." % \
                (i, states[i].temperature(), states[i].pressure(), str(self.be[i])))

            # Store the statistics of the equilibrium calculation
            iterations[i] = result.optimum.iterations
            seconds[i] = result.optimum.time

        # Extract the calculation statistics from the auxiliary ndarray members
        self.result.equilibrium.iterations.vector()[:] = iterations
        self.result.equilibrium.seconds.vector()[:] = seconds

        # Total time spent on performing equilibrium calculations
        self.result.seconds_equilibrium = time.time() - tbegin

    def transportKineticFluidSpecies(self, field, dt):
        pass

    def react(self, field, dt):
        pass

    def step(self, field, dt):

        if not self.initialized:
            self.initialize(field)

        tbegin = time.time()

        self.phi = self.porosity.phi()
        self.transportEquilibriumFluidSpecies(field, dt)
        self.transportKineticFluidSpecies(field, dt)
        self.react(field, dt)

        self.result.time = time.time() - tbegin
        self.result.succeeded = True
        return self.result

    def iterations(self):
        if self.equilibrium != []:
            self.output.vector()[:] = [x.iterations for x in self.equilibrium]
        else:
            self.output.assign(Constant(0.0))
        self.output.rename('IterationsPerDOF', 'IterationsPerDOF')
        return self.output

    def times(self):
        if self.equilibrium != []:
            self.output.vector()[:] = [x.time for x in self.equilibrium]
        else:
            self.output.assign(Constant(0.0))
        self.output.rename('SecondsPerDOF', 'SecondsPerDOF')
        return self.output


class ChemicalTransport:
    def __init__(self, mesh, system, **kwargs):
        # Create references to the parameters
        self.mesh = mesh
        self.system = system
        self.partition = Partition(system)
        self.mobility = Mobility(system)

        self.V = FunctionSpace(mesh, 'CG', 1)

        # Create a Function instance for outputting
        self.output = Function(self.V)

        self.u = Function(self.V)

        self.phi = Function(self.V)

        # Initialise the dictionary of phase velocities and diffusion coefficients
        velocity  = kwargs.get("velocity", Constant(0.0))
        diffusion = kwargs.get("diffusion", Constant(0.0))

        pore_velocity = velocity/self.phi
        pore_velocity_norm = sqrt(dot(pore_velocity, pore_velocity))
        pore_diffusion = diffusion*pore_velocity_norm

        # Initialise the transport solver
        self.transport = TransportSolver(self.u, velocity=pore_velocity, diffusion=pore_diffusion)

        # Initialise the equilibrium solver
        self.equilibrium = EquilibriumSolver(system)
        self.equilibrium.setPartition(self.partition)

        # Initialise the number of species and elements
        self.num_species  = system.numSpecies()
        self.num_elements = system.numElements()

        # Initialise the dof map
        self.dofmap = self.V.dofmap()

        # Initialise the coordinates of the dofs
        self.coordinates = self.dofmap.tabulate_all_coordinates(mesh)
        self.coordinates = self.coordinates.reshape((-1, mesh.geometry().dim()))

        # Initialise the number of degree-of-freedoms
        self.num_dofs = len(self.dofmap.dofs())

        # Initialise the chemical fluid of every degree-of-freedom
        self.states = [ChemicalState(system) for i in range(self.num_dofs)]

        # Initialise the arrays of molar amounts of each element for each degree of freedom
        self.be  = numpy.zeros((self.num_dofs, self.num_elements)) # in the equilibrium partition
        self.bef = numpy.zeros((self.num_dofs, self.num_elements)) # in the equilibrium-fluid partition
        self.bes = numpy.zeros((self.num_dofs, self.num_elements)) # in the equilibrium-solid partition

        # Flag that indicates if the solver has been initial

    def setPartition(self, partition):
        self.partition = partition
        self.equilibrium.setPartition(partition)

    def setMobility(self, mobility):
        self.mobility = mobility

    def setInitialCondition(self, state):
        for k in range(self.num_dofs):
            self.states[k].assign(state)
        self._updatePorosity()

    def setInitialConditionOnRegion(self, state, region):
        for k in range(self.num_dofs):
            if region(self.coordinates[k]):
                self.states[k].assign(state)
        self._updatePorosity()

    def setBoundaryCondition(self, state, boundary):
        self.bc = ChemicalDirichletBC(self.V, state, self.mobility, boundary)

    def setTemperatures(self, temperatures):
        for k in range(self.num_dofs):
            self.states[k].setTemperature(temperatures[k])

    def setPressures(self, pressures):
        for k in range(self.num_dofs):
            self.states[k].setPressure(pressures[k])

    def initialize(self):
        # Initialise the indices of the equilibrium and kinetic species
        self.ispecies_e  = self.partition.indicesEquilibriumSpecies()
        self.ispecies_k  = self.partition.indicesKineticSpecies()

        # Initialise the indices of the fluid and solid species
        self.ispecies_f  = self.mobility.indicesFluidSpecies()
        self.ispecies_s  = self.mobility.indicesSolidSpecies()

        # Initialise the indices of the equilibrium-fluid and equilibrium-solid species
        self.ispecies_ef = sorted(set(self.ispecies_e) & set(self.ispecies_f))
        self.ispecies_es = sorted(set(self.ispecies_e) & set(self.ispecies_s))

        # Initialise the indices of the kinetic-fluid and kinetic-solid species
        self.ispecies_kf = sorted(set(self.ispecies_k) & set(self.ispecies_f))
        self.ispecies_ks = sorted(set(self.ispecies_k) & set(self.ispecies_s))

    def _updatePorosity(self):
        porosities = numpy.empty(self.num_dofs)
        for k in range(self.num_dofs):
            T = self.states[k].temperature()
            P = self.states[k].pressure()
            n = self.states[k].speciesAmounts()
            v = self.system.phaseVolumes(T, P, n).val()
            isolid_phases = self.mobility.indicesSolidPhases()
            solid_volume = sum([v[i] for i in isolid_phases])
            porosities[k] = 1 - solid_volume
        self.phi.vector()[:] = porosities

    def transportEquilibriumFluidSpecies(self, dt):
        # Calculate the amounts of elements in the equilibrium-fluid partition
        elementAmountsInSpecies(self.states, self.ispecies_ef, self.bef)

        # Calculate the amounts of elements in the equilibrium-solid partition
        elementAmountsInSpecies(self.states, self.ispecies_es, self.bes)

        # Iterate over all elements in the equilibrium-fluid partition and tranport them
        for j in range(self.num_elements):
            # Create the boundary condition for the current element
            bc_element = self.bc.elementDirichletBC(j, self.ispecies_ef, self.phi)

            # Set initial condition for the transport equation of the current element
            self.u.vector()[:] = numpy.array(self.bef[:, j])

            # Transport the current element of the equilibrium-fluid partition
            self.transport.step(self.u, dt, bc_element)

            # Extract the result to the array of element amounts `bef`
            self.bef[:, j][:] = self.u.vector().get_local()

        # Ensure the new amounts of elements in the equilibrium-fluid partition is positive
        for k in range(self.num_dofs):
            for j in range(self.num_elements):
                self.bef[k, j] = max(1e-16, self.bef[k, j])

        # Calculate the new amounts of elements in the equilibrium partition
        numpy.add(self.bes, self.bef, self.be)

        # Equilibrate the chemical states for every degree-of-freedom
        self.equilibrate()

    def equilibrate(self):
        tbegin = time.time()
        self.equilibrium.solve(self.be, self.states)
        self.result.time_equilibrium = time.time() - tbegin

    def transportKineticFluidSpecies(self, dt):
        pass

    def react(self, dt):
        pass

    def step(self, dt):
        tbegin = time.time()

        self.result = ChemicalTransportResult()
        self.transportEquilibriumFluidSpecies(dt)
        self.transportKineticFluidSpecies(dt)
        self.react(dt)

        self.result.time = time.time() - tbegin
        self.result.succeeded = True
        return self.result

    def n(self, species):
        ispecies = self.system.indexSpeciesWithError(species)
        outvec = self.output.vector()
        amounts = numpy.empty(self.num_dofs)
        for k in range(self.num_dofs):
            amounts[k] = self.states[k].speciesAmount(ispecies)
        outvec[:] = amounts
        self.output.rename(species, species)
        return self.output

    def iterations(self):
        if self.equilibrium.result.iterations != []:
            self.output.vector()[:] = self.equilibrium.result.iterations
        else:
            self.output.assign(Constant(0.0))
        self.output.rename('IterationsPerDOF', 'IterationsPerDOF')
        return self.output

    def times(self):
        if self.equilibrium.result.time != []:
            self.output.vector()[:] = self.equilibrium.result.time
        else:
            self.output.assign(Constant(0.0))
        self.output.rename('SecondsPerDOF', 'SecondsPerDOF')
        return self.output

    def porosity(self):
        self.output.assign(self.phi)
        self.output.rename('Porosity', 'Porosity')
        return self.output

    def volume(self):
        volumes = numpy.empty(self.num_dofs)
        for k in range(self.num_dofs):
            T = self.states[k].temperature()
            P = self.states[k].pressure()
            n = self.states[k].speciesAmounts()
            v = self.system.phaseVolumes(T, P, n).val()
            volumes[k] = sum(v)
        self.output.vector()[:] = volumes
        self.output.rename('Volume', 'Volume')
        return self.output

    def ph(self):
        ln10 = math.log(10)
        iH = self.system.indexSpecies('H+')
        phvals = numpy.empty(self.num_dofs)
        for k in range(self.num_dofs):
            T = self.states[k].temperature()
            P = self.states[k].pressure()
            n = self.states[k].speciesAmounts()
            ln_a = self.system.lnActivities(T, P, n).val()
            phvals[k] = -ln_a[iH]/ln10
        self.output.vector()[:] = phvals
        self.output.rename('pH', 'pH')
        return self.output