/*
 * Reaktor is a C++ library for computational reaction modelling.
 *
 * Copyright (C) 2014 Allan Leal
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "ChemicalSystem.hpp"

// C++ includes
#include <algorithm>
#include <set>
#include <map>

// Reaktor includes
#include <Reaktor/Common/Exception.hpp>
#include <Reaktor/Core/Phase.hpp>
#include <Reaktor/Core/Species.hpp>
#include <Reaktor/Utils/SetUtils.hpp>

namespace Reaktor {
namespace internal {

auto inline unknownElementError(const std::string& element) -> void
{
    Exception exception;
    exception.error << "Cannot continue execution with an unknown element";
    exception.reason << "The element " << element << " was not specified in the chemical system.";
    raise(exception);
}

auto inline unknownSpeciesError(const std::string& species) -> void
{
    Exception exception;
    exception.error << "Cannot continue execution with an unknown species";
    exception.reason << "The species " << species << " was not specified in the chemical system.";
    raise(exception);
}

auto inline unknownPhaseError(const std::string& phase) -> void
{
    Exception exception;
    exception.error << "Cannot continue execution with an unknown phase";
    exception.reason << "The phase " << phase << " was not specified in the chemical system.";
    raise(exception);
}

} /* namespace internal */

using namespace internal;

class ChemicalSystem::Impl
{
private:
    /// The phases in the system
    std::vector<Phase> phases$;

    /// The chemical species in the system
    std::vector<Species> m_species;

    /// The chemical elements in the system
    std::vector<std::string> m_elements;

    /// The elemental formula matrix of the chemical species
    Matrix m_formula_matrix;

    /// The species map object that stores the indices of the species
    std::map<std::string, Index> m_species_map;

    /// The mapping from a chemical species to the index of the phase where it belongs
    Indices m_species_phase_map;

    /// The mapping from a phase to the indices of the chemical species that belongs to it
    std::vector<Indices> m_phase_species_map;

    /// The mapping from a chemical species to the indices of the chemical elements that compose it
    std::vector<Indices> m_species_element_map;

    /// The mapping from a chemical elements to the indices of the chemical species that contains it
    std::vector<Indices> m_element_species_map;

public:
    Impl()
    {}

    Impl(const std::vector<Phase>& phases)
    : phases$(phases)
    {
        initialiseSpecies();
        initialiseElements();
        initialiseFormulaMatrix();
        initialiseSpeciesMap();
        initialiseSpeciesPhaseMap();
        initialisePhaseSpeciesMap();
        initialiseSpeciesElementMap();
        initialiseElementSpeciesMap();
    }

    auto numPhases() const -> unsigned
    {
        return phases$.size();
    }

    auto numSpecies() const -> unsigned
    {
        return m_species.size();
    }

    auto numElements() const -> unsigned
    {
        return m_elements.size();
    }

    auto phase(const Index& idx) const -> const Phase&
    {
        return phases$[idx];
    }

    auto phases() const -> const std::vector<Phase>&
    {
        return phases$;
    }

    auto phasesNames() const -> std::vector<std::string>
    {
        std::vector<std::string> names;
        names.reserve(numPhases());
        for(const Phase& phase : phases$)
            names.push_back(phase.name());
        return names;
    }

    auto species(const Index& idx) const -> const Species&
    {
        return m_species[idx];
    }

    auto species(const std::string& name) const -> const Species&
    {
        return species(idxSpeciesWithError(name));
    }

    auto species() const -> const std::vector<Species>&
    {
        return m_species;
    }

    auto speciesNames() const -> std::vector<std::string>
    {
        return names(m_species);
    }

    auto speciesCharges() const -> Vector
    {
        return charges(m_species);
    }

    auto element(const Index& idx_element) const -> const std::string&
    {
        return m_elements[idx_element];
    }

    auto elements() const -> const std::vector<std::string>&
    {
        return m_elements;
    }

    auto formulaMatrix() const -> const Matrix&
    {
        return m_formula_matrix;
    }

    auto idxSpeciesInPhase(const Index& idx) const -> const Indices&
    {
        return m_phase_species_map[idx];
    }

    auto idxElementsInSpecies(const Index& idx) const -> const Indices&
    {
        return m_species_element_map[idx];
    }

    auto mapSpeciesToPhase() const -> const Indices&
    {
        return m_species_phase_map;
    }

    auto mapPhaseToSpecies() const -> const std::vector<Indices>&
    {
        return m_phase_species_map;
    }

    auto mapSpeciesToElements() const -> const std::vector<Indices>&
    {
        return m_species_element_map;
    }

    auto mapElementToSpecies() const -> const std::vector<Indices>&
    {
        return m_element_species_map;
    }

    auto idxElement(const std::string& name) const -> Index
    {
        return find(name, m_elements);
    }

    auto idxElements(const std::vector<std::string>& names) const -> Indices
    {
        return find(names, m_elements);
    }

    auto idxSpecies(const std::string& species) const -> Index
    {
        return m_species_map.find(species)->second;
    }

    auto idxSpecies(const std::vector<std::string>& names) const -> Indices
    {
        Indices indices;
        indices.reserve(names.size());
        for(const std::string name : names)
            indices.push_back(idxSpeciesWithError(name));
        return indices;
    }

    auto idxPhase(const std::string& name) const -> Index
    {
        auto same_name = [&](const Phase& phase) { return phase.name() == name; };

        return std::find_if(phases$.begin(), phases$.end(), same_name) - phases$.begin();
    }

    auto idxPhases(const std::vector<std::string>& names) const -> Indices
    {
        Indices indices;
        indices.reserve(names.size());
        for(const std::string name : names)
            indices.push_back(idxPhase(name));
        return indices;
    }

    auto idxPhaseWithSpecies(Index ispecies) const -> Index
    {
        return m_species_phase_map[ispecies];
    }

    auto idxPhaseWithSpecies(const std::string& species) const -> Index
    {
        const Index ispecies = idxSpeciesWithError(species);
        return idxPhaseWithSpecies(ispecies);
    }

    auto containsSpecies(const std::string& species) const -> bool
    {
        return idxSpecies(species) < numSpecies();
    }

    auto containsElement(const std::string& element) const -> bool
    {
        return idxElement(element) < numElements();
    }

    auto containsPhase(const std::string& phase) const -> bool
    {
        return idxPhase(phase) < numPhases();
    }

    auto idxElementWithError(const std::string& element) const -> Index
    {
        Index idxelement = idxElement(element);
        if(idxelement >= numElements())
            unknownElementError(element);
        return idxelement;
    }

    auto idxSpeciesWithError(const std::string& species) const -> Index
    {
        Index idxspecies = idxSpecies(species);
        if(idxspecies >= numSpecies())
            unknownSpeciesError(species);
        return idxspecies;
    }

    auto idxPhaseWithError(const std::string& phase) const -> Index
    {
        Index idxphase = idxPhase(phase);
        if(idxphase >= numPhases())
            unknownPhaseError(phase);
        return idxphase;
    }

    auto chemicalPotentials(double T, double P) const -> Vector
    {
        Vector mu0(numSpecies());
        unsigned offset = 0;
        for(const Phase& phase : phases$)
        {
            const unsigned Np = phase.numSpecies();
            mu0.segment(offset, Np) = phase.chemicalPotentials(T, P);
            offset += Np;
        }

        return mu0;
    }

    auto molarFractions(const Vector& n) const -> Vector
    {
        // The number of species in the system
        const auto N = numSpecies();

        // The molar fractions of the species
        Vector x = zeros(N);

        // An offset variable to determine how many species to process
        unsigned offset = 0;
        for(const Phase& phase : phases$)
        {
            // The number of species in the current phase
            const auto Np = phase.numSpecies();

            // The molar abundance of the species in the current phase
            const Vector np = n.segment(offset, Np);

            // Calculate the molar fractions in the current phase
            x.segment(offset, Np) = phase.molarFractions(np);

            // Increment the offset variable by the number of species in the current phase
            offset += Np;
        }

        return x;
    }

    auto concentrations(const Vector& n) const -> Vector
    {
        // The number of species in the system
        const auto N = numSpecies();

        // The concentrations of the species
        Vector c = zeros(N);

        // An offset variable to determine how many species to process
        unsigned offset = 0;
        for(const Phase& phase : phases$)
        {
            // The number of species in the current phase
            const auto Np = phase.numSpecies();

            // The molar abundance of the species in the current phase
            const Vector np = n.segment(offset, Np);

            // Calculate the concentrations in the current phase
            c.segment(offset, Np) = phase.concentrations(np);

            // Increment the offset variable by the number of species in the current phase
            offset += Np;
        }

        return c;
    }

    auto activities(double T, double P, const Vector& n) const -> PartialVector
    {
        const unsigned N = numSpecies();

        PartialVector res;
        func(res) = zeros(N);
        grad(res) = zeros(N, N);

        unsigned offset = 0;
        for(const Phase& phase : phases$)
        {
            const unsigned Np = phase.numSpecies();
            const Vector np = n.segment(offset, Np);
            auto pair = phase.activities(T, P, np);

            func(res).segment(offset, Np) = func(pair);
            grad(res).block(offset, offset, Np, Np) = grad(pair);

            offset += Np;
        }

        return res;
    }

    auto initialiseSpecies() -> void
    {
        // Collect the chemical species from all phases
        m_species.clear();
        for(const Phase& phase : phases$)
            m_species.insert(m_species.end(),
                phase.species().begin(), phase.species().end());
    }

    auto initialiseElements() -> void
    {
        // Collect the chemical elements from the species in the system
        std::set<std::string> set_elements;
        for(const Species& iter : m_species)
            for(auto pair : iter.elements())
                set_elements.insert(std::get<0>(pair));
        m_elements.assign(set_elements.begin(), set_elements.end());
    }

    auto initialiseFormulaMatrix() -> void
    {
        // Assemble the elemental formula matrix of the chemical species
        m_formula_matrix = zeros(m_elements.size(), m_species.size());
        for(unsigned i = 0; i < m_elements.size(); ++i)
            for(unsigned j = 0; j < m_species.size(); ++j)
                m_formula_matrix(i, j) = m_species[j].elementAtoms(m_elements[i]);
    }

    auto initialiseSpeciesMap() -> void
    {
        for(unsigned i = 0; i < numSpecies(); ++i)
            m_species_map[m_species[i].name()] = i;
    }

    auto initialiseSpeciesPhaseMap() -> void
    {
        // Initialise the mapping from the index of a species to the index of the phase where it belongs
        m_species_phase_map.resize(m_species.size());
        for(unsigned idx_phase = 0; idx_phase < phases$.size(); ++idx_phase)
            for(const std::string& species_name : phases$[idx_phase].species())
                m_species_phase_map[m_species_map[species_name]] = idx_phase;
    }

    auto initialisePhaseSpeciesMap() -> void
    {
        // Initialise the mapping from the index of a phase to the indices of the species that belongs to it
        for(const Phase& phase : phases$)
            m_phase_species_map.push_back(idxSpecies(names(phase.species())));
    }

    auto initialiseSpeciesElementMap() -> void
    {
        // Initialise the mapping from the index of a species to the indices of the elements that compose it
        m_species_element_map.resize(m_species.size());
        for(unsigned i = 0; i < m_species.size(); ++i)
            for(unsigned e = 0; e < m_elements.size(); ++e)
                if(m_formula_matrix(e, i) != 0.0)
                    m_species_element_map[i].push_back(e);
    }

    auto initialiseElementSpeciesMap() -> void
    {
        // Initialise the mapping from the index of a element to the indices of the species that contains it
        m_element_species_map.resize(m_elements.size());
        for(unsigned i = 0; i < m_species.size(); ++i)
            for(unsigned e = 0; e < m_elements.size(); ++e)
                if(m_formula_matrix(e, i) != 0.0)
                    m_element_species_map[e].push_back(i);
    }
};

ChemicalSystem::ChemicalSystem()
: pimpl(new Impl())
{}

ChemicalSystem::ChemicalSystem(const std::vector<Phase>& phases)
: pimpl(new Impl(phases))
{}

auto ChemicalSystem::numPhases() const -> unsigned
{
	return pimpl->numPhases();
}

auto ChemicalSystem::numSpecies() const -> unsigned
{
    return pimpl->numSpecies();
}

auto ChemicalSystem::numElements() const -> unsigned
{
	return pimpl->numElements();
}

auto ChemicalSystem::phase(const Index& idx_phase) const -> const Phase&
{
	return pimpl->phase(idx_phase);
}

auto ChemicalSystem::phases() const -> const std::vector<Phase>&
{
	return pimpl->phases();
}

auto ChemicalSystem::phasesNames() const -> std::vector<std::string>
{
    return pimpl->phasesNames();
}

auto ChemicalSystem::species(const Index& idx_species) const -> const Species&
{
	return pimpl->species(idx_species);
}

auto ChemicalSystem::species(const std::string& name) const -> const Species&
{
    return pimpl->species(name);
}

auto ChemicalSystem::species() const -> const std::vector<Species>&
{
	return pimpl->species();
}

auto ChemicalSystem::speciesNames() const -> std::vector<std::string>
{
    return pimpl->speciesNames();
}

auto ChemicalSystem::speciesCharges() const -> Vector
{
    return pimpl->speciesCharges();
}

auto ChemicalSystem::element(const Index& idx_element) const -> const std::string&
{
    return pimpl->element(idx_element);
}

auto ChemicalSystem::elements() const -> const std::vector<std::string>&
{
    return pimpl->elements();
}

auto ChemicalSystem::formulaMatrix() const -> const Matrix&
{
	return pimpl->formulaMatrix();
}

auto ChemicalSystem::idxSpeciesInPhase(const Index& idx_phase) const -> const Indices&
{
    return pimpl->idxSpeciesInPhase(idx_phase);
}

auto ChemicalSystem::idxElementsInSpecies(const Index& idx_species) const -> const Indices&
{
    return pimpl->idxElementsInSpecies(idx_species);
}

auto ChemicalSystem::mapSpeciesToPhase() const -> const Indices&
{
    return pimpl->mapSpeciesToPhase();
}

auto ChemicalSystem::mapPhaseToSpecies() const -> const std::vector<Indices>&
{
    return pimpl->mapPhaseToSpecies();
}

auto ChemicalSystem::mapSpeciesToElements() const -> const std::vector<Indices>&
{
    return pimpl->mapSpeciesToElements();
}

auto ChemicalSystem::mapElementToSpecies() const -> const std::vector<Indices>&
{
    return pimpl->mapElementToSpecies();
}

auto ChemicalSystem::idxElement(const std::string& element) const -> Index
{
    return pimpl->idxElement(element);
}

auto ChemicalSystem::idxElements(const std::vector<std::string>& elements) const -> Indices
{
    return pimpl->idxElements(elements);
}

auto ChemicalSystem::idxSpecies(const std::string& species) const -> Index
{
    return pimpl->idxSpecies(species);
}

auto ChemicalSystem::idxSpecies(const std::vector<std::string>& species) const -> Indices
{
    return pimpl->idxSpecies(species);
}

auto ChemicalSystem::idxPhase(const std::string& phase) const -> Index
{
	return pimpl->idxPhase(phase);
}

auto ChemicalSystem::idxPhases(const std::vector<std::string>& phases) const -> Indices
{
    return pimpl->idxPhases(phases);
}

auto ChemicalSystem::idxPhaseWithSpecies(Index ispecies) const -> Index
{
    return pimpl->idxPhaseWithSpecies(ispecies);
}

auto ChemicalSystem::idxPhaseWithSpecies(const std::string& species) const -> Index
{
	return pimpl->idxPhaseWithSpecies(species);
}

auto ChemicalSystem::containsSpecies(const std::string& species) const -> bool
{
    return pimpl->containsSpecies(species);
}

auto ChemicalSystem::containsElement(const std::string& element) const -> bool
{
    return pimpl->containsElement(element);
}

auto ChemicalSystem::containsPhase(const std::string& phase) const -> bool
{
    return pimpl->containsPhase(phase);
}

auto ChemicalSystem::idxElementWithError(const std::string& element) const -> Index
{
    return pimpl->idxElementWithError(element);
}

auto ChemicalSystem::idxSpeciesWithError(const std::string& species) const -> Index
{
    return pimpl->idxSpeciesWithError(species);
}

auto ChemicalSystem::idxPhaseWithError(const std::string& phase) const -> Index
{
    return pimpl->idxPhaseWithError(phase);
}

auto ChemicalSystem::chemicalPotentials(double T, double P) const -> Vector
{
    return pimpl->chemicalPotentials(T, P);
}

auto ChemicalSystem::molarFractions(const Vector& n) const -> Vector
{
    return pimpl->molarFractions(n);
}

auto ChemicalSystem::concentrations(const Vector& n) const -> Vector
{
    return pimpl->concentrations(n);
}

auto ChemicalSystem::activities(double T, double P, const Vector& n) const -> PartialVector
{
    return pimpl->activities(T, P, n);
}

auto operator<<(std::ostream& out, const ChemicalSystem& system) -> std::ostream&
{
    out << "Phases:" << std::endl;
    const auto& phases = system.phases();
    for(unsigned i = 0; i < phases.size(); ++i)
    {
        out << "  " << phases[i].name() << " Phase" << std::endl;
        out << "    ";
        const auto& species = phases[i].species();
        for(unsigned i = 0; i < species.size(); ++i)
        {
            const auto name = species[i].name();
            const auto idx  = system.idxSpecies(name);

            out << (i > 0 ? ", " : "") << idx << ":" << name;
        }
    }

	out << std::endl;

	out << "Species:" << std::endl;
	const auto& species = system.species();
	for(unsigned i = 0; i < species.size(); ++i)
	    out << (i > 0 ? ", " : "") << i << ":" << species[i].name();

	out << std::endl;

	out << "Elements:" << std::endl;
	const auto& elements = system.elements();
    for(unsigned i = 0; i < elements.size(); ++i)
        out << (i > 0 ? ", " : "") << i << ":" << elements[i];

    out << std::endl;

    out << "Formula Matrix:" << std::endl;
	out << system.formulaMatrix();

	return out;
}

} // namespace Reaktor
