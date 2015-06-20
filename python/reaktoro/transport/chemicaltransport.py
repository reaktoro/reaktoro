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
        assert type(velocity) is list, \
            'Expecting a list of velocity fields for each fluid phase.'
        assert len(velocity) == self.num_fluid_phases, \
            'There are %d fluid phases, but only %d velocity fields were given.' \
                % (self.num_fluid_phases, len(velocity))
        self.velocity = velocity
        self.initialized = False


    def setDiffusion(self, diffusion):
        assert type(diffusion) is list, \
            'Expecting a list of diffusion coefficients for each fluid phase.'
        assert len(diffusion) == self.num_fluid_phases, \
            'There are %d fluid phases, but only %d diffusion coefficients were given.' \
                % (self.num_fluid_phases, len(diffusion))
        self.diffusion = diffusion
        self.initialized = False


    def setSource(self, source):
        self.source = source
        self.initialized = False


    def setPartition(self, partition):
        if type(partition) is str:
            self.partition.set(partition)
        else:
            self.partition = partition


    def addBoundaryCondition(self, state, boundary):
        self.boundary_conditions.append((state, boundary))


    def setTemperatures(self, temperatures):
        if type(temperatures) in [int, float]:
            temperatures = [temperatures]*self.num_dofs
        assert hasattr(temperatures, '__len__'), 'Expecting a list of temperatures.'
        assert len(temperatures) == self.num_dofs, 'Expecting a list of \
            temperatures of length %d, but given list has length %d.' \
                % (self.num_dofs, len(temperatures))
        for i in range(self.num_dofs):
            self.states[i].setTemperature(temperatures[i])


    def setPressures(self, pressures):
        if type(pressures) in [int, float]:
            pressures = [pressures]*self.num_dofs
        assert hasattr(pressures, '__len__'), 'Expecting a list of pressures.'
        assert len(pressures) == self.num_dofs, 'Expecting a list of \
            pressures of length %d, but given list has length %d.' \
                % (self.num_dofs, len(pressures))
        for i in range(self.num_dofs):
            self.states[i].setPressure(pressures[i])


    def initialize(self, field):

        self.initialized = True

        # Check if the user provided any boundary conditions
        if self.boundary_conditions is []:
            RuntimeError('Failed to initialize ChemicalTransportSolver. \
                No boundary conditions have been provided.')

        # Initialize the function space used in the ChemicalField instance
        self.function_space = field.functionSpace()

        # Initialize the list of global `dof` indices in the current process
        self.dofs = self.function_space.dofmap().dofs()

        # Initialize the number of degrees-of-freedom in the current process
        self.num_dofs = len(self.dofs)

        # Initialize the ChemicalDirichletBC instances
        self.bcs = [ChemicalDirichletBC(self.function_space, state, self.mobility, boundary) \
                    for (state, boundary) in self.boundary_conditions]

        # The auxiliary Function instance used for transport steps
        self.u = Function(self.function_space)

        # The auxiliary Function instance used for outputting
        self.output = Function(self.function_space)

        # Initialize the chemical equilibrium solver
        self.equilibrium = EquilibriumSolver(system)
        self.equilibrium.setPartition(self.partition)

        # Initialize the indices of the equilibrium and kinetic species
        self.ispecies_e  = self.partition.indicesEquilibriumSpecies()
        self.ispecies_k  = self.partition.indicesKineticSpecies()

        # Initialize the indices of the fluid and solid species
        self.ispecies_f  = self.mobility.indicesFluidSpeciesInEachFluidPhase()
        self.ispecies_s  = self.mobility.indicesSolidSpecies()

        # Initialize the indices of the equilibrium-fluid and equilibrium-solid species
        self.ispecies_ef = [sorted(set(self.ispecies_e) & set(indices)) for indices in self.ispecies_f]
        self.ispecies_es = sorted(set(self.ispecies_e) & set(self.ispecies_s))

        # Initialize the indices of the kinetic-fluid and kinetic-solid species
        self.ispecies_kf = [sorted(set(self.ispecies_k) & set(indices)) for indices in self.ispecies_f]
        self.ispecies_ks = sorted(set(self.ispecies_k) & set(self.ispecies_s))

        # Initialize the indices of fluid and solid phases
        self.iphases_f = mobility.indicesFluidPhases()
        self.iphases_s = mobility.indicesSolidPhases()

        # Initialize the number of fluid and solid phases
        self.num_fluid_phases = len(self.iphases_f)
        self.num_solid_phases = len(self.iphases_s)

        # Initialize the arrays of element amounts for each fluid phase in
        # the equilibrium-fluid partition for each degree-of-freedom in the function space
        self.bef = [numpy.zeros((self.num_elements, self.num_dofs)) \
                        for i in xrange(self.num_fluid_phases)]

        # Initialize the array of element amounts in the equilibrium
        # partition for each degree-of-freedom in the function space
        self.be = numpy.zeros((self.num_elements, self.num_dofs))

        # Initialize the dolfin Function instances for the saturation field of each fluid phase
        self.saturation = [Function(self.function_space) for i in self.iphases_f]

        # Initialize the numpy array for the saturation field of each fluid phase
        self.saturation_array = numpy.zeros((self.num_fluid_phases, self.num_dofs))

        # Initialize the dolfin Function instance for the porosity field
        self.porosity = Function(self.function_space)

        # Initialize the numpy array for the porosity field
        self.porosity_array = numpy.zeros(self.num_dofs)

        # Define the dolfin forms for the pore velocity of each fluid phase
        self.pore_velocity = [self.velocity[i]/(self.porosity*self.saturation[i]) \
            for i in xrange(self.num_fluid_phases)]

        # Create transport solvers for each fluid phase
        self.transport = [TransportSolver() for i in xrange(self.num_fluid_phases)]

        # Set the pore velocities and diffusion coefficients for each transport solver
        for i in xrange(self.num_fluid_phases):
            self.transport[i].setVelocity(self.pore_velocity[i])
            self.transport[i].setDiffusion(self.diffusion[i])

        # Initialize the ChemicalTransportResult instance
        self.result = ChemicalTransportResult()
        self.result.equilibrium.iterations = Function(self.function_space)
        self.result.equilibrium.seconds = Function(self.function_space)
        self.result.kinetics.timesteps = Function(self.function_space)
        self.result.kinetics.seconds = Function(self.function_space)

        self.result.equilibrium._iterations = numpy.empty(self.num_dofs)
        self.result.equilibrium._seconds = numpy.empty(self.num_dofs)
        self.result.kinetics._timesteps = numpy.empty(self.num_dofs)
        self.result.kinetics._seconds = numpy.empty(self.num_dofs)

        self.result.equilibrium.iterations.rename('EquilibriumIterationsPerDOF', 'EquilibriumIterationsPerDOF')
        self.result.equilibrium.seconds.rename('EquilibriumSecondsPerDOF', 'EquilibriumSecondsPerDOF')
        self.result.kinetics.timesteps.rename('KineticsTimeStepsPerDOF', 'KineticsTimeStepsPerDOF')
        self.result.kinetics.seconds.rename('KineticsSecondsPerDOF', 'KineticsSecondsPerDOF')


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


    def elementDirichletBC(self, ielement, iphase):
        ispecies = self.ispecies_ef[iphase] # the indices of equilibrium species in phase `iphase`
        return [bc.elementDirichletBC(ielement, ispecies, self.porosity) for bc in self.bcs]


    def speciesDirichletBC(self, ispecies):
        return [bc.speciesDirichletBC(ispecies, self.porosity) for bc in self.bcs]


    def transportEquilibriumElementInFluidPhase(self, ielement, iphase, field, dt):
        # The Dirichlet boundary conditions element `ielement` in fluid phase `iphase`
        bcs = self.elementDirichletBC(ielement, iphase)

        # Set the boundary conditions in the tranport solver w.r.t. fluid phase `iphase`
        self.transport[iphase].setBoundaryConditions(bcs)

        # Calculate the element amounts in fluid phase `iphase` of the equilibrium-fluid partition
        field.elementAmountsInSpecies(self.ispecies_ef[iphase], self.bef[iphase])

        # Set initial condition for the transport equation of the current element
        self.u.vector()[:] = self.bef[iphase][ielement]

        # Transport the current element of the equilibrium-fluid partition
        self.transport[iphase].step(self.u, dt)

        # Extract the result to the array of element amounts `bef`
        self.bef[iphase][ielement] = self.u.vector().get_local()


    def transportEquilibriumElementsInFluidPhases(self, field, dt):

        # Iterate over all elements in all fluid phases in the equilibrium-fluid
        # partition and tranport them
        for iphase in xrange(self.num_fluid_phases):
            for ielement in xrange(self.num_elements):
                self.transportEquilibriumElementInFluidPhase(ielement, iphase, field, dt)


    def transportKineticSpeciesInFluidPhases(self, field, dt):
        pass


    def equilibrate(self, field):
        # Start timing the equilibrate step
        tbegin = time.time()

        # Define auxiliary bindings
        states = field.states()
        iterations = self.result.equilibrium._iterations
        seconds = self.result.equilibrium._seconds

        # Calculate the element amounts in the equilibrium partition.
        # Compute the contribution from the equilibrium-solid partition
        field.elementAmountsInSpecies(self.ispecies_es, self.be)

        # Compute the contribution from the equilibrium-fluid partition
        for bef in self.bef:
            self.be += bef

        # Compute the equilibrium state in the equilibrium
        for i in xrange(self.num_dofs):
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

    def react(self, field, dt):
        pass

    def step(self, field, dt):
        # Check if the chemical transport solver has been initialized
        if not self.initialized:
            self.initialize(field)

        # Start timing the transport step calculation
        tbegin = time.time()

        # Update the porosity field and saturation fields of each fluid phase
        self.updatePorositySaturation(field)

        # Transport the elements in the equilibrium-fluid partition
        self.transportEquilibriumElementsInFluidPhases(field, dt)

        # Transport the species in the kinetic-fluid partition
        self.transportKineticSpeciesInFluidPhases(field, dt)

        # Equilibrate the equilibrium-solid and equilibrium-fluid partitions
        self.equilibrate(field)

        # React the equilibrium and kinetic partitions over a time of `dt`
        self.react(field, dt)

        # Compute the time elapsed for the chemical transport step
        self.result.time = time.time() - tbegin
