# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright (C) 2014-2020 Allan Leal
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library. If not, see <http://www.gnu.org/licenses/>.

from reaktoro import *
from pytest import approx, raises
from numpy import array


def names(objects):
    """Extract the names of each object in a list."""
    return [o.name() for o in objects]


def range_species_in_phase(system, iphase):
    """Return the range of species indices in a phase."""
    return range(system.indexFirstSpeciesInPhase(iphase),
        system.indexFirstSpeciesInPhase(iphase) + system.numSpeciesInPhase(iphase))


def test_chemical_system():
    """Test function for class ChemicalSystem."""

    editor = ChemicalEditor()
    editor.addAqueousPhase("H2O(l) H+ OH- HCO3- CO2(aq) CO3--".split())
    editor.addGaseousPhase("H2O(g) CO2(g)".split())
    editor.addMineralPhase("Graphite")
    system = ChemicalSystem(editor)

    # The number of elements, species, and phases
    Ne = system.numElements()
    Ns = system.numSpecies()
    Np = system.numPhases()

    # A sensible value for temperature (in K)
    T = 300

    # A sensible value for pressure (in Pa)
    P = 1e5

    # A sensible array of species amounts
    n = array([55, 1e-7, 1e-7, 0.1, 0.5, 0.01, 1.0, 0.001, 1.0])

    # The formula matrix of the system
    A = system.formulaMatrix()

    # The amounts of the elements in the system
    b = A @ n

    # The formula matrices of each phase in the system
    Ap = [ A[:, range_species_in_phase(system, j)] for j in range(Np) ]

    # The amounts of the species of each phase in the system
    np = [ n[range_species_in_phase(system, j)] for j in range(Np) ]

    # The amounts of the elements of each phase in the system
    bp = [ Ap[j] @ np[j] for j in range(Np) ]

    # Calculate the chemical properties of the state
    properties = system.properties(T, P, n)

    # -------------------------------------------------------------------------
    # Check methods ChemicalSystem::num(Elements|Species|Phases)
    # -------------------------------------------------------------------------
    assert Ne == len(system.elements())
    assert Ns == len(system.species())
    assert Np == len(system.phases())

    # -------------------------------------------------------------------------
    # Check method ChemicalSystem::numSpeciesInPhase
    # -------------------------------------------------------------------------
    for j in range(Np):
        assert system.numSpeciesInPhase(j) == len(system.phase(j).species())

    # -------------------------------------------------------------------------
    # Check method ChemicalSystem::formulaMatrix
    # -------------------------------------------------------------------------
    assert system.formulaMatrix() == approx(array([
        [ s.elementCoefficient(e.name()) for s in system.species()]
            for e in system.elements() ]))

    # -------------------------------------------------------------------------
    # Check method ChemicalSystem::element
    # -------------------------------------------------------------------------
    for i, element in enumerate(system.elements()):
        assert system.element(i) == element

    for i, element in enumerate(system.elements()):
        assert system.element(element.name()) == element

    # -------------------------------------------------------------------------
    # Check method ChemicalSystem::species
    # -------------------------------------------------------------------------
    for i, species in enumerate(system.species()):
        assert system.species(i) == species

    # Check method ChemicalSystem::species(name)
    for species in system.species():
        assert system.species(species.name()) == species

    # -------------------------------------------------------------------------
    # Check method ChemicalSystem::phase
    # -------------------------------------------------------------------------
    for i, phase in enumerate(system.phases()):
        assert system.phase(i) == phase

    for i, phase in enumerate(system.phases()):
        assert system.phase(phase.name()) == phase

    # -------------------------------------------------------------------------
    # Check method ChemicalSystem::indexElement
    # -------------------------------------------------------------------------
    for i, element in enumerate(system.elements()):
        assert system.indexElement(element.name()) == i

    for i, element in enumerate(system.elements()):
        assert system.indexElementWithError(element.name()) == i

    with raises(RuntimeError):
        assert system.indexElementWithError("Aa")

    # -------------------------------------------------------------------------
    # Check method ChemicalSystem::indexSpecies
    # -------------------------------------------------------------------------
    for i, species in enumerate(system.species()):
        assert system.indexSpecies(species.name()) == i

    for i, species in enumerate(system.species()):
        assert system.indexSpeciesWithError(species.name()) == i

    with raises(RuntimeError):
        assert system.indexSpeciesWithError("AaBb2")

    # -------------------------------------------------------------------------
    # Check method ChemicalSystem::indexSpeciesAny
    # -------------------------------------------------------------------------
    for i in range(Ns):
        name = system.species(i).name()
        assert system.indexSpeciesAny(["AaBb", name]) == system.indexSpecies(name)

    with raises(RuntimeError):
        assert system.indexSpeciesAnyWithError(["AaBb2", "AaCc3"])

    # -------------------------------------------------------------------------
    # Check method ChemicalSystem::indexPhase
    # -------------------------------------------------------------------------
    for i, phase in enumerate(system.phases()):
        assert system.indexPhase(phase.name()) == i

    for i, phase in enumerate(system.phases()):
        assert system.indexPhaseWithError(phase.name()) == i

    with raises(RuntimeError):
        assert system.indexPhaseWithError("WrongName")

    # -------------------------------------------------------------------------
    # Check method ChemicalSystem::indexPhaseWithSpecies
    # -------------------------------------------------------------------------
    for iphase, phase in enumerate(system.phases()):
        for ispecies, species in enumerate(phase.species()):
            i = system.indexSpecies(species.name())
            assert system.indexPhaseWithSpecies(i) == iphase

    # -------------------------------------------------------------------------
    # Check method ChemicalSystem::indexFirstSpeciesInPhase
    # -------------------------------------------------------------------------
    for iphase, phase in enumerate(system.phases()):
        i = system.indexSpecies(phase.species(0).name())
        assert system.indexFirstSpeciesInPhase(iphase) == i

    # -------------------------------------------------------------------------
    # Check method ChemicalSystem::indicesElements
    # -------------------------------------------------------------------------
    assert system.indicesElements(names(system.elements())) == approx(range(Ne))

    # -------------------------------------------------------------------------
    # Check method ChemicalSystem::indicesElementsInSpecies
    # -------------------------------------------------------------------------
    for ispecies, species in enumerate(system.species()):
        ielements = [system.indexElement(e.name()) for e in species.elements()]
        assert system.indicesElementsInSpecies(ispecies) == approx(ielements)

    # -------------------------------------------------------------------------
    # Check method ChemicalSystem::indicesSpecies
    # -------------------------------------------------------------------------
    assert system.indicesSpecies(names(system.species())) == approx(range(Ns))

    # -------------------------------------------------------------------------
    # Check method ChemicalSystem::indicesSpeciesInPhases
    # -------------------------------------------------------------------------
    assert system.indicesSpeciesInPhases(range(Np)) == approx(range(Ns))

    for iphase, phase in enumerate(system.phases()):
        assert system.indicesSpeciesInPhases([iphase]) == \
            system.indicesSpecies(names(phase.species()))

    # -------------------------------------------------------------------------
    # Check method ChemicalSystem::indicesPhases
    # -------------------------------------------------------------------------
    assert system.indicesPhases(names(system.phases())) == approx(range(Np))

    # -------------------------------------------------------------------------
    # Check method ChemicalSystem::indicesPhasesWithSpecies
    # -------------------------------------------------------------------------
    assert system.indicesPhasesWithSpecies(array(range(Ns))) == approx(range(Np))

    # -------------------------------------------------------------------------
    # Check method ChemicalSystem::indicesFluidPhases
    # -------------------------------------------------------------------------
    assert system.indicesFluidPhases() == \
        approx([i for i, phase in enumerate(system.phases()) if phase.isFluid()])

    # -------------------------------------------------------------------------
    # Check method ChemicalSystem::indicesFluidSpecies
    # -------------------------------------------------------------------------
    assert system.indicesFluidSpecies() == \
        approx(system.indicesSpeciesInPhases(system.indicesFluidPhases()))

    # -------------------------------------------------------------------------
    # Check method ChemicalSystem::indicesSolidPhases
    # -------------------------------------------------------------------------
    assert system.indicesSolidPhases() == \
        approx([i for i, phase in enumerate(system.phases()) if phase.isSolid()])

    # -------------------------------------------------------------------------
    # Check method ChemicalSystem::indicesSolidSpecies
    # -------------------------------------------------------------------------
    assert system.indicesSolidSpecies() == \
        approx(system.indicesSpeciesInPhases(system.indicesSolidPhases()))

    # -------------------------------------------------------------------------
    # Check method ChemicalSystem::elementAmounts
    # -------------------------------------------------------------------------
    assert system.elementAmounts(n) == approx(b)

    # -------------------------------------------------------------------------
    # Check method ChemicalSystem::elementAmountsInPhase
    # -------------------------------------------------------------------------
    for j in range(Np):
        assert system.elementAmountsInPhase(j, n) == approx(bp[j])

    # -------------------------------------------------------------------------
    # Check method ChemicalSystem::elementAmountsInSpecies
    # -------------------------------------------------------------------------
    assert system.elementAmountsInSpecies([], n) == approx(0.0)
    assert system.elementAmountsInSpecies(range(Ns), n) == approx(b)

    # -------------------------------------------------------------------------
    # Check method ChemicalSystem::elementAmount
    # -------------------------------------------------------------------------
    for j in range(Ne):
        assert system.elementAmount(j, n) == approx(b[j])

    # -------------------------------------------------------------------------
    # Check method ChemicalSystem::elementAmountInPhase
    # -------------------------------------------------------------------------
    for ielement, iphase in ((x, y) for x in range(Ne) for y in range(Np)):
        assert system.elementAmountInPhase(ielement, iphase, n) == approx(bp[iphase][ielement])

    # -------------------------------------------------------------------------
    # Check method ChemicalSystem::elementAmountInSpecies
    # -------------------------------------------------------------------------
    for ielement in range(Ne):
        assert system.elementAmountInSpecies(ielement, [], n) == approx(0.0)
        assert system.elementAmountInSpecies(ielement, range(Ns), n) == approx(b[ielement])

    # Check the usage system.properties(T, P, n).someProperty() works
    assert all(system.properties(T, P, n).phaseVolumes().val == properties.phaseVolumes().val)
    assert all(system.properties(T, P, n).lnActivities().val == properties.lnActivities().val)
    assert all(system.properties(T, P, n).chemicalPotentials().val == properties.chemicalPotentials().val)
