import math
import numpy as npy
from dolfin import *
from reaktor import *
from geoflow.transport import Transport
from geoflow.equilibrium import Equilibrium
from geoflow.mobility import Mobility
from numpy import empty
import time

def elementAmounts(states, result):
    for k in range(len(states)):
        result[k] = states[k].elementAmounts()

def elementAmountsInSpecies(states, ispecies_list, result):
    for k in range(len(states)):
        bspecies = states[k].elementAmountsInSpecies(ispecies_list)
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
    def __init__(self, V, state, boundary):
        self.state = state
        self.boundary = boundary
        self.values = Function(V)
        self.dirichlet = DirichletBC(V, 1, boundary)
        self.dirichlet.apply(self.values.vector())
        self.dofs = npy.where(self.values.vector() == 1)[0]
        self.dirichlet = DirichletBC(V, self.values, boundary)

    def elementDirichletBC(self, ielement, ispecies_list, porosity):
        bval = self.state.elementAmountInSpecies(ielement, ispecies_list)
        phi = porosity.vector()[self.dofs]
        self.values.vector()[self.dofs] = phi*bval
        return self.dirichlet

    def speciesDirichletBC(self, ispecies, porosity):
        nval = self.state.speciesAmount(ispecies)
        phi = porosity.vector()[self.dofs]
        self.values.vector()[self.dofs] = phi*nval
        return self.dirichlet


class ChemicalTransportResult:
    def __init__(self):
        # The flag that indicates if the chemical transport calculations succeeded
        self.succeeded = False

        # The total time of the chemical transport calculation
        self.time = 0

        # The total time spent on equilibrium calculations
        self.time_equilibrium = 0

        # The total time spent on finite element assembly operations
        self.time_assemble = 0

        # The total time spent on linear systems
        self.time_linear_systems = 0

    def __iadd__(self, other):
        self.succeeded            = self.succeeded and other.succeeded
        self.time                += other.time
        self.time_equilibrium    += other.time_equilibrium
        self.time_assemble       += other.time_assemble
        self.time_linear_systems += other.time_linear_systems
        return self

class ChemicalTransport:
    def __init__(self, mesh, system, partition, mobility, **kwargs):
        # Create references to the parameters
        self.mesh = mesh
        self.system = system
        self.partition = partition
        self.mobility = mobility

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
        self.transport = Transport(self.u, velocity=pore_velocity, diffusion=pore_diffusion)

        # Initialise the equilibrium solver
        self.equilibrium = Equilibrium(system, partition, **kwargs)

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

        # Initialise the number of species and elements
        self.num_species    = system.numSpecies()
        self.num_elements   = system.numElements()

        # Initialise the indices of the equilibrium and kinetic species
        self.num_species_e  = len(self.ispecies_e)
        self.num_species_k  = len(self.ispecies_k)

        # Initialise the indices of the fluid and solid species
        self.num_species_f  = len(self.ispecies_f)
        self.num_species_s  = len(self.ispecies_s)

        # Initialise the indices of the equilibrium-fluid and equilibrium-solid species
        self.num_species_ef = len(self.ispecies_ef)
        self.num_species_es = len(self.ispecies_es)

        # Initialise the indices of the kinetic-fluid and kinetic-solid species
        self.num_species_kf = len(self.ispecies_kf)
        self.num_species_ks = len(self.ispecies_ks)

        # Initialise the dof map
        self.dofmap = self.V.dofmap()

        # Initialise the coordinates of the dofs
        self.coordinates = self.dofmap.tabulate_all_coordinates(mesh)
        self.coordinates = self.coordinates.reshape((-1, mesh.geometry().dim()))

        # Initialise the number of degree-of-freedoms
        self.num_dofs = len(self.dofmap.dofs())

        # Initialise the chemical state of every degree-of-freedom
        self.states = [ChemicalState(system) for i in range(self.num_dofs)]

        # Initialise the arrays of molar amounts of each element for each degree of freedom
        self.be  = npy.zeros((self.num_dofs, self.num_elements)) # in the equilibrium partition
        self.bef = npy.zeros((self.num_dofs, self.num_elements)) # in the equilibrium-fluid partition
        self.bes = npy.zeros((self.num_dofs, self.num_elements)) # in the equilibrium-solid partition

        # Flags to indicate if properties have been updated
        self.updated_velocity  = True
        self.updated_diffusion = True
        self.updated_source    = True

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
        fluid = self.mobility.fluid(state)
        self.bc = ChemicalDirichletBC(self.V, fluid, boundary)

    def setTemperatures(self, temperatures):
        for k in range(self.num_dofs):
            self.states[k].setTemperature(temperatures[k])

    def setPressures(self, pressures):
        for k in range(self.num_dofs):
            self.states[k].setPressure(pressures[k])

    def n(self, species):
        ispecies = self.system.indexSpeciesWithError(species)
        outvec = self.output.vector()
        amounts = npy.empty(self.num_dofs)
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

#     def porosity(self):
#         porosities = empty(self.num_dofs)
#         for k in range(self.num_dofs):
#             T = self.states[k].temperature()
#             P = self.states[k].pressure()
#             n = self.states[k].speciesAmounts()
#             v = self.system.phaseVolumes(T, P, n).val()
#             fluid_volume = sum([v[i] for i in self.mobility.indicesFluidPhases()])
#             total_volume = sum(v)
#             porosities[k] = fluid_volume/total_volume
#         self.output.vector()[:] = porosities
#         self.output.rename('Porosity', 'Porosity')
#         return self.output

    def porosity(self):
        self.output.assign(self.phi)
        self.output.rename('Porosity', 'Porosity')
        return self.output

    def volume(self):
        volumes = npy.empty(self.num_dofs)
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
        phvals = npy.empty(self.num_dofs)
        for k in range(self.num_dofs):
            T = self.states[k].temperature()
            P = self.states[k].pressure()
            n = self.states[k].speciesAmounts()
            ln_a = self.system.lnActivities(T, P, n).val()
            phvals[k] = -ln_a[iH]/ln10
        self.output.vector()[:] = phvals
        self.output.rename('pH', 'pH')
        return self.output


    def _updatePorosity(self):
        porosities = npy.empty(self.num_dofs)
        for k in range(self.num_dofs):
            T = self.states[k].temperature()
            P = self.states[k].pressure()
            n = self.states[k].speciesAmounts()
            v = self.system.phaseVolumes(T, P, n).val()
            isolid_phases = self.mobility.indicesSolidPhases()
            solid_volume = sum([v[i] for i in isolid_phases])
            porosities[k] = 1 - solid_volume
        self.phi.vector()[:] = porosities

    def _transportEquilibriumFluidSpecies(self, dt):
        # Calculate the amounts of elements in the equilibrium-fluid partition
        elementAmountsInSpecies(self.states, self.ispecies_ef, self.bef)

        # Calculate the amounts of elements in the equilibrium-solid partition
        elementAmountsInSpecies(self.states, self.ispecies_es, self.bes)

        # Iterate over all elements in the equilibrium-fluid partition and tranport them
        for j in range(self.num_elements):
            # Create the boundary condition for the current element
            bc_element = self.bc.elementDirichletBC(j, self.ispecies_ef, self.phi)

            # Set initial condition for the transport equation of the current element
            self.u.vector()[:] = npy.array(self.bef[:, j])

            # Transport the current element of the equilibrium-fluid partition
            self.transport.step(self.u, dt, bc_element)

            # Extract the result to the array of element amounts `bef`
            self.bef[:, j][:] = self.u.vector().get_local()

        # Ensure the new amounts of elements in the equilibrium-fluid partition is positive
        for k in range(self.num_dofs):
            for j in range(self.num_elements):
                self.bef[k, j] = max(1e-16, self.bef[k, j])

        # Calculate the new amounts of elements in the equilibrium partition
        npy.add(self.bes, self.bef, self.be)

        # Equilibrate the chemical states for every degree-of-freedom
        self._equilibrate()

    def _equilibrate(self):
        tbegin = time.time()
        self.equilibrium.solve(self.be, self.states)
        self.result.time_equilibrium = time.time() - tbegin

    def _transportKineticFluidSpecies(self, dt):
        pass

    def _react(self, dt):
        pass

    def step(self, dt):
        tbegin = time.time()

        self.result = ChemicalTransportResult()
        self._transportEquilibriumFluidSpecies(dt)
        self._transportKineticFluidSpecies(dt)
        self._react(dt)

        self.result.time = time.time() - tbegin
        self.result.succeeded = True
        return self.result
