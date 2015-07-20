// Reaktoro is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "Interface.hpp"

// C++ includes
#include <map>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Element.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Core/Species.hpp>

namespace Reaktoro {
namespace {

/// Return the Element instances and their stoichiometries that compose a species
auto elementsInSpecies(const Interface& interface, std::vector<Element> elements, Index ispecies) -> std::map<Element, double>
{
    std::map<Element, double> res;
    for(unsigned i = 0; i < interface.numElements(); ++i)
        if(interface.elementStoichiometry(ispecies, i) != 0.0)
            res.emplace(elements[i], interface.elementStoichiometry(ispecies, i));
    return res;
}

/// Return the Species instances that compose a phase
auto speciesInPhase(const Interface& interface, std::vector<Species> species, Index iphase) -> std::vector<Species>
{
    const unsigned ifirst = interface.indexFirstSpeciesInPhase(iphase);
    const unsigned nspecies = interface.numSpeciesInPhase(iphase);
    const auto begin = species.begin() + ifirst;
    const auto end = begin + nspecies;
    return std::vector<Species>(begin, end);
}

/// Convert a string to a PhaseReferenceState value
auto convertReferenceState(std::string ref) -> PhaseReferenceState
{
    if(ref == "IdealGas") return PhaseReferenceState::IdealGas;
    if(ref == "IdealSolution") return PhaseReferenceState::IdealSolution;
    RuntimeError("Could not set the standard reference state of the phase.",
        "The given string `" + ref + "` is neither IdealGas nor IdealSolution.");
}

/// Return a PhaseThermoModel function for a phase
auto phaseThermoModel(Interface* interface, Index iphase) -> PhaseThermoModel
{
    const unsigned ifirst = interface->indexFirstSpeciesInPhase(iphase);
    const unsigned nspecies = interface->numSpeciesInPhase(iphase);

    double T_old = 0.0, P_old = 0.0;

    PhaseThermoModelResult res(nspecies);

    PhaseThermoModel f = [=](double T, double P) mutable
    {
        if(T != T_old || P != P_old)
        {
            T_old = T;
            P_old = P;
            interface->set(T, P);
        }

        res.standard_partial_molar_gibbs_energies.val     = rows(interface->standardMolarGibbsEnergies(), ifirst, nspecies);
        res.standard_partial_molar_enthalpies.val         = rows(interface->standardMolarEnthalpies(), ifirst, nspecies);
        res.standard_partial_molar_volumes.val            = rows(interface->standardMolarVolumes(), ifirst, nspecies);
        res.standard_partial_molar_heat_capacities_cp.val = rows(interface->standardMolarHeatCapacitiesConstP(), ifirst, nspecies);
        res.standard_partial_molar_heat_capacities_cv.val = rows(interface->standardMolarHeatCapacitiesConstV(), ifirst, nspecies);
        return res;
    };

    return f;
}

/// Return a PhaseChemicalModel function for a phase
auto phaseChemicalModel(Interface* interface, Index iphase) -> PhaseChemicalModel
{
    const unsigned ifirst = interface->indexFirstSpeciesInPhase(iphase);
    const unsigned nspecies = interface->numSpeciesInPhase(iphase);

    double T_old = 0.0, P_old = 0.0;
    Vector n_old;

    PhaseChemicalModelResult res(nspecies);

    PhaseChemicalModel f = [=](double T, double P, const Vector& n) mutable
    {
        if(T != T_old || P != P_old || n != n_old)
        {
            T_old = T;
            P_old = P;
            n_old = n;
            interface->set(T, P, n);
        }

        res.ln_activity_coefficients.val        = rows(interface->lnActivityCoefficients(), ifirst, nspecies);
        res.ln_activities.val                   = rows(interface->lnActivities(), ifirst, nspecies);
        res.molar_volume.val                    = interface->phaseMolarVolumes()[iphase];
        res.residual_molar_gibbs_energy.val     = interface->phaseResidualMolarGibbsEnergies()[iphase];
        res.residual_molar_enthalpy.val         = interface->phaseResidualMolarEnthalpies()[iphase];
        res.residual_molar_heat_capacity_cp.val = interface->phaseResidualMolarHeatCapacitiesConstP()[iphase];
        res.residual_molar_heat_capacity_cv.val = interface->phaseResidualMolarHeatCapacitiesConstV()[iphase];
        return res;
    };

    return f;
}

} // namespace

auto Interface::phaseMolarVolumes() const -> Vector
{
    return zeros(numPhases());
}

auto Interface::phaseResidualMolarGibbsEnergies() const -> Vector
{
    return zeros(numPhases());
}

auto Interface::phaseResidualMolarEnthalpies() const -> Vector
{
    return zeros(numPhases());
}

auto Interface::phaseResidualMolarHeatCapacitiesConstP() const -> Vector
{
    return zeros(numPhases());
}

auto Interface::phaseResidualMolarHeatCapacitiesConstV() const -> Vector
{
    return zeros(numPhases());
}

auto Interface::formulaMatrix() const -> Matrix
{
    const unsigned E = numElements();
    const unsigned N = numSpecies();
    Matrix A(E, N);
    for(unsigned i = 0; i < N; ++i)
        for(unsigned j = 0; j < E; ++j)
            A(j, i) = elementStoichiometry(i, j);
    return A;
}

auto Interface::indexElement(std::string element) const -> Index
{
    const Index size = numElements();
    for(unsigned i = 0; i < size; ++i)
        if(elementName(i) == element)
            return i;
    return size;
}

auto Interface::indexSpecies(std::string species) const -> Index
{
    const Index size = numSpecies();
    for(unsigned i = 0; i < size; ++i)
        if(speciesName(i) == species)
            return i;
    return size;
}

auto Interface::indexPhase(std::string phase) const -> Index
{
    const Index size = numPhases();
    for(unsigned i = 0; i < size; ++i)
        if(phaseName(i) == phase)
            return i;
    return size;
}

auto Interface::indexPhaseWithSpecies(Index ispecies) const -> Index
{
    Assert(ispecies < numSpecies(), "Cannot get the index of the phase with a given species.",
        "The given species index `" + std::to_string(ispecies) + "` is out of range.");
    const Index size = numPhases();
    unsigned counter = 0;
    for(unsigned i = 0; i < size; ++i)
    {
        counter += numSpeciesInPhase(i);
        if(counter > ispecies) return i;
    }
    return size;
}

auto Interface::indexFirstSpeciesInPhase(Index iphase) const -> Index
{
    Assert(iphase < numPhases(), "Cannot get the index of first species in a given phase.",
        "The given phase index `" + std::to_string(iphase) + "` is out of range.");
    Index counter = 0;
    for(unsigned i = 0; i < iphase; ++i)
        counter += numSpeciesInPhase(iphase);
    return counter;
}

Interface::operator ChemicalSystem()
{
    const unsigned nelements = numElements();
    const unsigned nspecies = numSpecies();
    const unsigned nphases = numPhases();

    const Interface& interface = *this;

    // Create the Element instances
    std::vector<Element> elements(nelements);
    for(unsigned i = 0; i < nelements; ++i)
    {
        elements[i].setName(elementName(i));
        elements[i].setMolarMass(elementMolarMass(i));
    }

    // Create the Species instances
    std::vector<Species> species(nspecies);
    for(unsigned i = 0; i < nspecies; ++i)
    {
        species[i].setName(speciesName(i));
        species[i].setFormula(speciesName(i));
        species[i].setElements(elementsInSpecies(interface, elements, i));
    }

    // Create the Phase instances
    std::vector<Phase> phases(nphases);
    for(unsigned i = 0; i < nphases; ++i)
    {
        phases[i].setName(phaseName(i));
        phases[i].setSpecies(speciesInPhase(interface, species, i));
        phases[i].setReferenceState(convertReferenceState(phaseReferenceState(i)));
        phases[i].setThermoModel(phaseThermoModel(this, i));
        phases[i].setChemicalModel(phaseChemicalModel(this, i));
    }

    // Create the ChemicalSystem instance
    ChemicalSystem system(phases);

    return system;
}

Interface::operator ChemicalState()
{
    ChemicalSystem system = *this;
    ChemicalState state(system);
    state.setTemperature(temperature());
    state.setPressure(pressure());
    state.setSpeciesAmounts(speciesAmounts());
    return state;
}

} // namespace Reaktoro
