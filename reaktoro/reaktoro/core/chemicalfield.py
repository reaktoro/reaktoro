import numpy
from reaktoro.core import ChemicalState
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


    def _set_fill(self, state):
        for k in range(self.num_dofs):
            self.states[k].assign(state)

    def _set_where(self, state, where):
        for k in range(self.num_dofs):
            if where(self.coordinates[k]) is True:
                self.states[k].assign(state)

    def set(self, state, **kwargs):
        where = kwargs.get('where')
        if where is not None:
            self._set_where(state, where)
        else:
            self._set_fill(state)

    def setTemperatures(self, temperatures):
        for k in range(self.num_dofs):
            self.states[k].setTemperature(temperatures[k])

    def setPressures(self, pressures):
        for k in range(self.num_dofs):
            self.states[k].setPressure(pressures[k])

    def elementAmounts(self, b):
        for i in range(self.num_dofs):
            vec = self.states[i].elementAmounts()
            b[:, i] = vec.array()

    def elementAmountsInSpecies(self, indices, b):
        for i in range(self.num_dofs):
            vec = self.states[i].elementAmountsInSpecies(indices)
            b[:, i] = vec.array()

    def porosity(self):
        for k in range(self.num_dofs):
            state = self.states[k]
            T = state.temperature()
            P = state.pressure()
            n = state.speciesAmounts()
            v = system.phaseVolumes(T, P, n).val()
            fluid_volume = sum([v[i] for i in self.iphases_fluid])
            solid_volume = sum([v[i] for i in self.iphases_solid])
            self.values[k] = 1.0 - solid_volume
        self.porosity.vector()[:] = self.values
        return self.porosity


class ChemicalField(object):

    def __init__(self, system, function_space):
        self.pimpl = _ChemicalField(system, function_space)

    def set(self, state, **kwargs):
        self.pimpl.set(state, kwargs)

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

    def functionSpace(self):
        return self.pimpl.function_space

