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

#include "Database.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringList.hpp>
#include <Reaktoro/Core/AggregateState.hpp>
#include <Reaktoro/Core/Element.hpp>
#include <Reaktoro/Core/Species.hpp>

namespace Reaktoro {
namespace {

/// Auxiliary type to store species with same aggregate state.
using SpeciesSameAggregateStateMap = std::unordered_map<AggregateState, std::vector<Species>>;

/// Return the values in a map in a vector container.
template<typename Map>
auto vectorize(const Map& map)
{
    using T = typename Map::mapped_type;
    std::vector<T> values;
    values.reserve(map.size());
    for(const auto& [_, value] : map)
        values.push_back(value);
    return values;
}

} // namespace

struct Database::Impl
{
    /// The map of species name to Species objects in the database.
    SpeciesMap species_map;

    /// The map of element name to Element objects in the database collected from the species.
    ElementMap element_map;

    /// The map from aggregate state to species with that aggregate state.
    SpeciesSameAggregateStateMap species_with_same_aggregate_state;

    /// The additional data in the database whose type is known at runtime only.
    std::any attached_data;

    /// Add a species in the database.
    auto addSpecies(Species species) -> void
    {
        // Ensure a unique name is used when storing this species
        auto name = species.name();
        auto unique_name = name;
        auto already_exists = true;
        while(already_exists) {
            already_exists = species_map.find(unique_name) != species_map.end();
            if(already_exists)
                unique_name = unique_name + "!"; // keep adding symbol ! to the name such as H2O, H2O!, H2O!!, H2O!!! if many H2O are given
        }

        // Replace name if not unique
        if(name != unique_name) {
            species = species.withName(unique_name);
            warning(true, "Species in a Database object should have unique names, but species ", name, " violates this rule. "
                "I'm changing its name to ", unique_name, " to fix this issue for you.");
        }

        species_map[species.name()] = std::move(species);
        for(auto&& [element, coeff] : species.elements())
            element_map[element.symbol()] = element;
        species_with_same_aggregate_state[species.aggregateState()].push_back(species);
    }

    /// Set additional data for the database whose type is known at runtime only.
    auto attachData(const std::any& data) -> void
    {
        attached_data = data;
    }

    /// Return all elements in the database.
    auto elements() const -> std::vector<Element>
    {
        return vectorize(element_map);
    }

    /// Return all species in the database.
    auto species() const -> std::vector<Species>
    {
        return vectorize(species_map);
    }

    /// Return the additional data in the database whose type is known at runtime only.
    auto attachedData() const -> const std::any&
    {
        return attached_data;
    }

    /// Return an element with given name if it exists in the database.
    auto elementWithName(std::string name) const -> std::optional<Element>
    {
        const auto iter = element_map.find(name);
        const auto found = iter != element_map.end();
        return found ? iter->second : std::optional<Element>{};
    }

    /// Return an element with given name if it exists in the database.
    auto speciesWithName(std::string name) const -> std::optional<Species>
    {
        const auto iter = species_map.find(name);
        const auto found = iter != species_map.end();
        return found ? iter->second : std::optional<Species>{};
    }

    /// Return all species in the database with given aggregate state.
    auto allSpeciesWithAggregateState(AggregateState state) const -> std::vector<Species>
    {
        const auto iter = species_with_same_aggregate_state.find(state);
        const auto found = iter != species_with_same_aggregate_state.end();
        return found ? iter->second : std::vector<Species>{};
    }

    /// Return all species in the database with given elements.
    /// @param symbols The symbols of the elements.
    auto allSpeciesWithElements(const StringList& symbols) const -> std::vector<Species>
    {
        auto with_elements = [&](const Species& species)
        {
            for(auto [element, _] : species.elements())
                if(!contains(symbols, element.name()))
                    return false;
            return true;
        };

        std::vector<Species> result;
        for(const auto& [_, species] : species_map)
            if(with_elements(species))
                result.push_back(species);
        return result;
    }

    /// Return true if the database contains an element with given name.
    /// @param name The name of the element.
    auto containsElement(std::string name) const -> bool
    {
        return element_map.find(name) != element_map.end();
    }

    /// Return true if the database contains a species with given name.
    /// @param name The name of the species.
    auto containsSpecies(std::string name) const -> bool
    {
        return species_map.find(name) != species_map.end();
    }
};

Database::Database()
: pimpl(new Impl())
{}

auto Database::addSpecies(const Species& species) -> void
{
    pimpl->addSpecies(species);
}

auto Database::attachData(const std::any& data) -> void
{
    pimpl->attachData(data);
}

auto Database::elements() const -> std::vector<Element>
{
    return pimpl->elements();
}

auto Database::species() const -> std::vector<Species>
{
    return pimpl->species();
}

auto Database::attachedData() const -> const std::any&
{
    return pimpl->attachedData();
}

auto Database::elementWithName(std::string name) const -> std::optional<Element>
{
    return pimpl->elementWithName(name);
}

auto Database::speciesWithName(std::string name) const -> std::optional<Species>
{
    return pimpl->speciesWithName(name);
}

auto Database::allSpeciesWithAggregateState(AggregateState state) const -> std::vector<Species>
{
    return pimpl->allSpeciesWithAggregateState(state);
}

auto Database::allSpeciesWithElements(const StringList& symbols) const -> std::vector<Species>
{
    return pimpl->allSpeciesWithElements(symbols);
}

} // namespace Reaktoro

