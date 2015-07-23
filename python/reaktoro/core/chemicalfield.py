import numpy
from core import ChemicalState
from dolfin import *

class _ChemicalField(object):

    def __init__(self, system, mobility, function_space):
        # Initialize the ChemicalSystem instance
        self.system = system

        # Initialize the Mobility instance
        self.mobility = mobility

        # Initialize the function space where the chemical field is defined
        self.function_space = function_space

        # Initialize the mesh member
        self.mesh = function_space.mesh()

        # Initialize the dof map
        self.dofmap = self.function_space.dofmap()

        # Initialize the number of degree-of-freedoms
        self.num_dofs = len(self.dofmap.dofs())

        # Initialize the coordinates of the dofs
        self.coordinates = self.dofmap.tabulate_all_coordinates(self.mesh)
        self.coordinates = self.coordinates.reshape((-1, self.mesh.geometry().dim()))

        # Initialize the chemical state of every degree-of-freedom
        self.states = [ChemicalState(system) for i in range(self.num_dofs)]

        # Initialize the indices of fluid and solid phases
        self.iphases_fluid = mobility.indicesFluidPhases()
        self.iphases_solid = mobility.indicesSolidPhases()

        # Initialize the Function instances for the saturation field of each fluid phase
        self.s = [Function(function_space) for i in self.iphases_fluid]

        # Initialize the Function instance for the porosity field
        self.phi = Function(function_space)

        # Initialize the auxliary Function instance for output purposes
        self.out = Function(function_space)

        # Initialize the auxliary array for output purposes
        self.values = numpy.zeros(self.num_dofs)


    def fill(self, state):
        for k in xrange(self.num_dofs):
            self.states[k].assign(state)

    def set(self, state, region):
        for k in xrange(self.num_dofs):
            if region(self.coordinates[k]) is True:
                self.states[k].assign(state)

    def setTemperatures(self, temperatures):
        for k in xrange(self.num_dofs):
            self.states[k].setTemperature(temperatures[k])

    def setPressures(self, pressures):
        for k in xrange(self.num_dofs):
            self.states[k].setPressure(pressures[k])

    def elementAmounts(self, b):
        for i in xrange(self.num_dofs):
            vec = self.states[i].elementAmounts()
            b[:, i] = vec.array()

    def elementAmountsInSpecies(self, indices, b):
        for i in xrange(self.num_dofs):
            vec = self.states[i].elementAmountsInSpecies(indices)
            b[:, i] = vec.array()

    def porosity(self):
        for k in xrange(self.num_dofs):
            state = self.states[k]
            T = state.temperature()
            P = state.pressure()
            n = state.speciesAmounts()
            properties = self.system.properties(T, P, n)
            v = properties.phaseVolumes().val
            solid_volume = sum([v[i] for i in self.iphases_solid])
            self.values[k] = 1.0 - solid_volume
        self.out.vector()[:] = self.values
        self.out.rename('Porosity', 'Porosity')
        return self.out

    def n(self, species):
        ispecies = self.system.indexSpecies(species)
        for k in xrange(self.num_dofs):
            self.values[k] = self.states[k].speciesAmount(ispecies)
        self.out.vector()[:] = self.values
        self.out.rename(species, species)
        return self.out

    def elementAmount(self, element, phase):
        ielement = self.system.indexElement(element)
        iphase = self.system.indexPhase(phase)
        for k in xrange(self.num_dofs):
            self.values[k] = self.states[k].elementAmountInPhase(ielement, iphase)
        self.out.vector()[:] = self.values
        self.out.rename(element, element)
        return self.out

#     def ph(self, species):
#         ispecies = self.system.indexSpecies(species)
#         for k in xrange(self.num_dofs):
#             self.values[k] = self.states[k].speciesAmount(ispecies)
#         self.out.vector()[:] = self.values
#         return self.out


class ChemicalField(object):

    def __init__(self, system, mobility, function_space):
        self.pimpl = _ChemicalField(system, mobility, function_space)

    def fill(self, state):
        self.pimpl.fill(state)

    def set(self, state, region):
        self.pimpl.set(state, region)

    def setTemperatures(self, temperatures):
        self.pimpl.setTemperatures(temperatures)

    def setPressures(self, pressures):
        self.pimpl.setPressures(pressures)

    def elementAmounts(self, b):
        self.pimpl.elementAmounts(b)

    def elementAmountsInSpecies(self, indices, b):
        self.pimpl.elementAmountsInSpecies(indices, b)

    def states(self):
        return self.pimpl.states

    def system(self):
        return self.pimpl.system

    def mobility(self):
        return self.pimpl.mobility

    def functionSpace(self):
        return self.pimpl.function_space

    def porosity(self):
        return self.pimpl.porosity()

    def n(self, species):
        return self.pimpl.n(species)

    def elementAmount(self, element, phase):
        return self.pimpl.elementAmount(element, phase)

