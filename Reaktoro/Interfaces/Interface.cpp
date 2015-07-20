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

    PhaseThermoModel f = [=](double T, double P) mutable
    {
        interface->set(T, P);
        PhaseThermoModelResult res(nspecies);
        for(unsigned i = 0; i < nspecies; ++i)
        {
            const unsigned j = ifirst + i;
            res.standard_partial_molar_gibbs_energies[i] = interface->standardMolarGibbsEnergy(j);
            res.standard_partial_molar_enthalpies[i] = interface->standardMolarEnthalpy(j);
            res.standard_partial_molar_volumes[i] = interface->standardMolarVolume(j);
            res.standard_partial_molar_heat_capacities_cp[i] = interface->standardMolarHeatCapacityConstP(j);
            res.standard_partial_molar_heat_capacities_cv[i] = interface->standardMolarHeatCapacityConstV(j);
        }
        return res;
    };

    return f;
}

/// Return a PhaseChemicalModel function for a phase
auto phaseChemicalModel(Interface* interface, Index iphase) -> PhaseChemicalModel
{
    const unsigned nspecies = interface->numSpecies();

    PhaseChemicalModel f = [=](double T, double P, const Vector& n)
    {
        interface->set(T, P, n);
        PhaseChemicalModelResult res(nspecies);
        res.ln_activity_coefficients.val = interface->lnActivityCoefficients(iphase);
        res.ln_activities.val = interface->lnActivities(iphase);
        res.molar_volume.val = interface->phaseMolarVolume(iphase);
        res.residual_molar_gibbs_energy.val = interface->phaseResidualMolarGibbsEnergy(iphase);
        res.residual_molar_enthalpy.val = interface->phaseResidualMolarEnthalpy(iphase);
        res.residual_molar_heat_capacity_cp.val = interface->phaseResidualMolarHeatCapacityConstP(iphase);
        res.residual_molar_heat_capacity_cv.val = interface->phaseResidualMolarHeatCapacityConstV(iphase);
        return res;
    };

    return f;
}

} // namespace

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
