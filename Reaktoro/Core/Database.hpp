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

#pragma once

// C++ includes
#include <any>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace Reaktoro {

// Forward declarations
class Element;
class Species;

/// Auxiliary type alias for an unordered map of element name to Element object.
using ElementMap = std::unordered_map<std::string, Element>;

/// Auxiliary type alias for an unordered map of species name to Species object.
using SpeciesMap = std::unordered_map<std::string, Species>;

/// The class used to store and retrieve data of chemical species.
/// @see Element, Species
/// @ingroup Databases
class Database
{
public:
    /// Construct a default Database instance.
    Database();

    /// Set all element in the database.
    auto setElements(const ElementMap& element_map) -> void;

    /// Set all species in the database.
    auto setSpecies(const SpeciesMap& species_map) -> void;

    /// Set additional data for the database whose type is known at runtime only.
    auto setData(const std::any& data) -> void;

    /// Add an element in the database.
    auto addElement(const Element& element) -> void;

    /// Add a species in the database.
    auto addSpecies(const Species& species) -> void;

    /// Return all elements in the database.
    auto elements() const -> std::vector<Element>;

    /// Return all species in the database.
    auto species() const -> std::vector<Species>;

    /// Return the additional data in the database whose type is known at runtime only.
    auto data() const -> const std::any&;

    /// Return an element with given name if it exists in the database.
    auto elementWithName(std::string name) const -> std::optional<Element>;

    /// Return an element with given name if it exists in the database.
    auto speciesWithName(std::string name) const -> std::optional<Species>;

    /// Return all species in the database with given type.
    auto speciesWithType(std::string type) const -> std::vector<Species>;

    /// Return all species in the database with given elements.
    /// @param elements The names of the elements.
    auto speciesWithElements(std::vector<std::string> elements) const -> std::vector<Species>;

    /// Return true if the database contains an element with given name.
    /// @param name The name of the element.
    auto containsElement(std::string name) const -> bool;

    /// Return true if the database contains a species with given name.
    /// @param name The name of the species.
    auto containsSpecies(std::string name) const -> bool;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

} // namespace Reaktoro
