import numpy
from dolfin import *

class _Saturation(object):

    def __init__(self, field, mobility):
        self.field = field
        self.mobility = mobility

        V = field.functionSpace()
        ifluid_phases = mobility.indicesFluidPhases()

        self.saturation = [Function(V) for i in ifluid_phases]

        self.dofmap = V.dofmap()
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


class Saturation(object):

    def __init__(self, field, mobility):
        self.pimpl = _Saturation(field, mobility)

    def phi(self):
        return self.pimpl.phi()
