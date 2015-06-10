

class Mobility:
    def __init__(self, system):
        self.system = system
        self._ifluid_phases  = []
        self._isolid_phases  = range(system.numPhases())
        self._ifluid_species = []
        self._isolid_species = range(system.numSpecies())

    def setFluidPhases(self, fluid_phases):
        phase_to_species = self.system.connectivity().phase_to_species
        for phase_name in fluid_phases:
            iphase = self.system.indexPhaseWithError(phase_name)
            self._ifluid_phases.append(iphase)
            self._ifluid_species += phase_to_species[iphase]
        prange = range(self.system.numPhases())
        srange = range(self.system.numSpecies())
        self._isolid_phases = set(prange) - set(self._ifluid_phases)
        self._isolid_species = set(srange) - set(self._ifluid_species)

    def setSolidPhases(self, solid_phases):
        fluid_phases = [phase.name() for phase in self.system.phases()]
        fluid_phases = set(fluid_phases) - set(solid_phases)
        self.setFluidPhases(fluid_phases)

    def indicesFluidPhases(self):
        return self._ifluid_phases

    def indicesSolidPhases(self):
        return self._isolid_phases

    def indicesFluidSpecies(self):
        return self._ifluid_species

    def indicesSolidSpecies(self):
        return self._isolid_species

    def fluid(self, state):
        res = state.clone()
        for i in self.indicesSolidPhases():
            res.setPhaseVolume(i, 0.0)
        res.setVolume(1.0)
        return res

    def solid(self, state):
        res = state.clone()
        for i in self.indicesFluidPhases():
            res.setPhaseVolume(i, 0.0)
        res.setVolume(1.0)
        return res

