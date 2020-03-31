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

// Forward declarations (class)
class Element;
class Species;
class StringList;

// Forward declarations (enum class)
enum class AggregateState;

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

    /// Add a species in the database.
    auto addSpecies(const Species& species) -> void;

    /// Attach data to this database whose type is known at runtime only.
    auto attachData(const std::any& data) -> void;

    /// Return all elements in the database.
    auto elements() const -> std::vector<Element>;

    /// Return all species in the database.
    auto species() const -> std::vector<Species>;

    /// Return the attached data to this database whose type is known at runtime only.
    auto attachedData() const -> const std::any&;

    /// Return an element with given name if it exists in the database.
    auto elementWithName(std::string name) const -> std::optional<Element>;

    /// Return a species with given name if it exists in the database.
    auto speciesWithName(std::string name) const -> std::optional<Species>;

    /// Return all species in the database with given aggregate state.
    auto allSpeciesWithAggregateState(AggregateState state) const -> std::vector<Species>;

    /// Return all species in the database with given elements.
    /// @param symbols The symbols of the elements.
    auto allSpeciesWithElements(const StringList& symbols) const -> std::vector<Species>;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

} // namespace Reaktoro
