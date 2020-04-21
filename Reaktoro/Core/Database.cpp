// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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

#include "Database.hpp"

// C++ includes
#include <deque>
#include <set>
#include <unordered_map>

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>

namespace Reaktoro {

struct Database::Impl
{
    /// The Species objects in the database.
    std::deque<Species> species;

    /// The Element objects in the database.
    std::set<Element> elements;

    /// The additional data in the database whose type is known at runtime only.
    std::any attached_data;

    /// The names of all species already in the database.
    std::set<std::string> species_names;

    /// The species in the database grouped in terms of their aggregate state
    std::unordered_map<AggregateState, std::deque<Species>> species_with_aggregate_state;

    /// Add a species in the database.
    auto addSpecies(Species newspecies) -> void
    {
        // Ensure a unique name is used when storing this species!
        auto name = newspecies.name();
        auto unique_name = name;
        auto already_exists = true;
        while(already_exists) {
            already_exists = species_names.find(unique_name) != species_names.end();
            if(already_exists)
                unique_name = unique_name + "!"; // keep adding symbol ! to the name such as H2O, H2O!, H2O!!, H2O!!! if many H2O are given
        }

        // Replace name if not unique.
        if(name != unique_name) {
            newspecies = newspecies.withName(unique_name);
            warning(true, "Species should have unique names in Database, but species ", name, " "
                "violates this rule. The unique name ", unique_name, " has been assigned instead.");
        }

        // Replace aggregate state to Aqueous if undefined. If this default
        // aggregate state option is not appropriate for a given scenario,
        // ensure then that the correct aggregate state of the species has
        // already been set in the Species object using method
        // Species::withAggregateState. This permits species names such as H2O,
        // CO2, H2, O2 as well as ions H+, OH-, HCO3- to be considered as
        // aqueous species, without requiring explicit aggregate state suffix
        // aq as in H2O(aq), CO2(aq), H2(aq). For all other species, ensure an
        // aggregate state suffix is provided, such as H2O(l), H2O(g),
        // CaCO3(s), Fe3O4(s).
        if(newspecies.aggregateState() == AggregateState::Undefined)
            newspecies = newspecies.withAggregateState(AggregateState::Aqueous);

        // Append the new Species object in the species container
        species.push_back(newspecies);

        // Update the list of unique species names.
        species_names.insert(newspecies.name());

        // Update the container of elements with the Element objects in this new species.
        for(auto&& [element, coeff] : newspecies.elements())
            elements.insert(element);

        // Add the new species in the group of species with same aggregate state
        species_with_aggregate_state[newspecies.aggregateState()].push_back(newspecies);
    }
};

Database::Database()
: pimpl(new Impl())
{}

Database::Database(const Database& other)
: pimpl(new Impl(*other.pimpl))
{}

Database::~Database()
{}

auto Database::operator=(Database other) -> Database&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto Database::addSpecies(const Species& species) -> void
{
    pimpl->addSpecies(species);
}

auto Database::attachData(const std::any& data) -> void
{
    pimpl->attached_data = data;
}

auto Database::elements() const -> ElementList
{
    const auto begin = pimpl->elements.begin();
    const auto end = pimpl->elements.end();
    return ElementList(begin, end);
}

auto Database::species() const -> SpeciesList
{
    const auto begin = pimpl->species.begin();
    const auto end = pimpl->species.end();
    return SpeciesList(begin, end);
}

auto Database::speciesWithAggregateState(AggregateState option) const -> SpeciesList
{
    auto it = pimpl->species_with_aggregate_state.find(option);
    if(it == pimpl->species_with_aggregate_state.end())
        return {};
    const auto begin = it->second.begin();
    const auto end = it->second.end();
    return SpeciesList(begin, end);
}

auto Database::attachedData() const -> const std::any&
{
    return pimpl->attached_data;
}

} // namespace Reaktoro

