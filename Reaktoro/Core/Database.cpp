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
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Core/Element.hpp>
#include <Reaktoro/Core/Species.hpp>

namespace Reaktoro {
namespace {

/// Auxiliary type alias for an unordered map of species category to SpeciesMap object.
using SpeciesGroupMap = std::unordered_map<std::string, SpeciesMap>;

/// Return values in a map in a vector.
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
    /// The map of element name to Element objects in the database.
    ElementMap element_map;

    /// The map of species name to Species objects in the database.
    SpeciesMap species_map;

    /// The map of species type to a a SpeciesMap object containing species with that type.
    SpeciesGroupMap species_group_with_type;

    /// Set all element in the database.
    auto setElements(const ElementMap& map) -> void
    {
        element_map = map;
    }

    /// Set all species in the database.
    auto setSpecies(const SpeciesMap& map) -> void
    {
        species_map = map;
        species_group_with_type.clear();
        for(const auto& [_, species] : species_map)
            species_group_with_type[species.type()][species.name()] = species;
    }

    /// Add an element in the database.
    auto addElement(const Element& element) -> void
    {
        element_map[element.name()] = element;
    }

    /// Add a species in the database.
    auto addSpecies(const Species& species) -> void
    {
        species_map[species.name()] = species;
        species_group_with_type[species.type()][species.name()] = species;
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

    /// Return all species in the database with given type.
    auto speciesWithType(std::string type) const -> std::vector<Species>
    {
        const auto iter = species_group_with_type.find(type);
        const auto found = iter != species_group_with_type.end();
        return found ? vectorize(iter->second) : std::vector<Species>{};
    }

    /// Return all species in the database with given elements.
    /// @param elements The names of the elements.
    auto speciesWithElements(std::vector<std::string> elements) const -> std::vector<Species>
    {
        auto with_elements = [&](const Species& species)
        {
            for(auto [element, _] : species.elements())
                if(!contained(element.name(), elements))
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

auto Database::setElements(const ElementMap& element_map) -> void
{
    pimpl->setElements(element_map);
}

auto Database::setSpecies(const SpeciesMap& species_map) -> void
{
    pimpl->setSpecies(species_map);
}

auto Database::addElement(const Element& element) -> void
{
    pimpl->addElement(element);
}

auto Database::addSpecies(const Species& species) -> void
{
    pimpl->addSpecies(species);
}

auto Database::elements() const -> std::vector<Element>
{
    return pimpl->elements();
}

auto Database::species() const -> std::vector<Species>
{
    return pimpl->species();
}

auto Database::elementWithName(std::string name) const -> std::optional<Element>
{
    return pimpl->elementWithName(name);
}

auto Database::speciesWithName(std::string name) const -> std::optional<Species>
{
    return pimpl->speciesWithName(name);
}

auto Database::speciesWithType(std::string type) const -> std::vector<Species>
{
    return pimpl->speciesWithType(type);
}

auto Database::speciesWithElements(std::vector<std::string> elements) const -> std::vector<Species>
{
    return pimpl->speciesWithElements(elements);
}

auto Database::containsElement(std::string name) const -> bool
{
    return pimpl->containsElement(name);
}

auto Database::containsSpecies(std::string name) const -> bool
{
    return pimpl->containsSpecies(name);
}

} // namespace Reaktoro

