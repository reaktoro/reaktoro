import numpy
from reaktoro.core import ChemicalState
from dolfin import *

class _ChemicalField(object):

    def __init__(self, system, function_space):
        # Initialise the chemical system
        self.system = system

        # Initialise the function space where the chemical field is defined
        self.function_space = function_space

        # Initialise the mesh member
        self.mesh = function_space.mesh()

        # Initialise the dof map
        self.dofmap = self.function_space.dofmap()

        # Initialise the number of degree-of-freedoms
        self.num_dofs = len(self.dofmap.dofs())

        # Initialise the coordinates of the dofs
        self.coordinates = self.dofmap.tabulate_all_coordinates(self.mesh)
        self.coordinates = self.coordinates.reshape((-1, self.mesh.geometry().dim()))

        # Initialise the chemical state of every degree-of-freedom
        self.states = [ChemicalState(system) for i in range(self.num_dofs)]

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
        for k in range(self.num_dofs):
            vec = self.states[k].elementAmounts()
            b[k] = vec.array()

    def elementAmountsInSpecies(self, indices, b):
        for k in range(self.num_dofs):
            vec = self.states[k].elementAmountsInSpecies(indices)
            b[k] = vec.array()


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

