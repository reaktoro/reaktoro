// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#include "Connectivity.hpp"

// Reaktoro includes
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>

namespace Reaktoro {

struct Connectivity::Impl
{
    /// The mapping from the index of an element to the indices of the species that contains it.
    std::vector<Indices> element_to_species;

    /// The mapping from the index of a species to the indices of the elements that it contains.
    std::vector<Indices> species_to_elements;

    /// The mapping from the index of a species to the index of the phase that contains it.
    Indices species_to_phase;

    /// The mapping from the index of a phase to the indices of the species that it contains.
    std::vector<Indices> phase_to_species;

    /// The mapping from the index of an element to the indices of the phases that contains it.
    std::vector<Indices> element_to_phases;

    /// The mapping from the index of a phase to the indices of the elements that it contains.
    std::vector<Indices> phase_to_elements;

    Impl()
    {}

    Impl(const ChemicalSystem& system)
    {
        const auto& elements = system.elements();
        const auto& species = system.species();
        const auto& phases = system.phases();

        const unsigned num_elements = elements.size();
        const unsigned num_species = species.size();
        const unsigned num_phases = phases.size();

        element_to_species.resize(num_elements);
        species_to_elements.resize(num_species);
        for(unsigned j = 0; j < num_elements; ++j)
            for(unsigned i = 0; i < num_species; ++i)
                if(species[i].elements().count(elements[j])) {
                    element_to_species[j].push_back(i);
                    species_to_elements[i].push_back(j); }

        species_to_phase.resize(num_species);
        phase_to_species.resize(num_phases);
        for(unsigned k = 0; k < num_phases; ++k)
            for(unsigned i = 0; i < num_species; ++i)
                if(contained(species[i], phases[k].species())) {
                    species_to_phase[i] = k;
                    phase_to_species[k].push_back(i);
                }

        element_to_phases.resize(num_elements);
        phase_to_elements.resize(num_phases);
        for(unsigned k = 0; k < num_phases; ++k)
            for(unsigned j = 0; j < num_elements; ++j)
                if(contained(elements[j], phases[k].elements())) {
                    element_to_phases[j].push_back(k);
                    phase_to_elements[k].push_back(j);
                }
    }
};

Connectivity::Connectivity()
: pimpl(new Impl())
{}

Connectivity::Connectivity(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

Connectivity::Connectivity(const Connectivity& other)
: pimpl(new Impl(*other.pimpl))
{}

Connectivity::~Connectivity()
{}

auto Connectivity::operator=(Connectivity other) -> Connectivity&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto Connectivity::indicesElementsInSpecies(Index ispecies) const -> const Indices&
{
    return pimpl->species_to_elements[ispecies];
}

auto Connectivity::indicesElementsInPhase(Index iphase) const -> const Indices&
{
    return pimpl->phase_to_elements[iphase];
}

auto Connectivity::indicesSpeciesInPhase(Index iphase) const -> const Indices&
{
    return pimpl->phase_to_species[iphase];
}

auto Connectivity::indicesSpeciesWithElement(Index ielement) const -> const Indices&
{
    return pimpl->element_to_species[ielement];
}

auto Connectivity::indicesPhasesWithElement(Index ielement) const -> const Indices&
{
    return pimpl->element_to_phases[ielement];
}

auto Connectivity::indexPhaseWithSpecies(Index ispecies) const -> Index
{
    return pimpl->species_to_phase[ispecies];
}

} // namespace Reaktoro
