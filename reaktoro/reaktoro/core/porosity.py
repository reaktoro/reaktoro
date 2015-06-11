import numpy
from dolfin import *

class _Porosity(object):

    def __init__(self, field, mobility):
        self.field = field
        self.mobility = mobility
        self.porosity = Function(field.functionSpace())
        self.dofmap = field.functionSpace().dofmap()
        self.num_dofs = len(self.dofmap.dofs())
        self.values = numpy.empty(self.num_dofs)

    def phi(self):
        states = self.field.states()
        system = self.field.system()
        isolid_phases = self.mobility.indicesSolidPhases()
        for k in range(self.num_dofs):
            T = states[k].temperature()
            P = states[k].pressure()
            n = states[k].speciesAmounts()
            v = system.phaseVolumes(T, P, n).val()
            solid_volume = sum([v[i] for i in isolid_phases])
            self.values[k] = 1.0 - solid_volume
        self.porosity.vector()[:] = self.values
        return self.porosity


class Porosity(object):

    def __init__(self, field, mobility):
        self.pimpl = _Porosity(field, mobility)

    def phi(self):
        return self.pimpl.phi()
