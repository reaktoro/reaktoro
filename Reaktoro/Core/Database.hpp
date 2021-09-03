// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/ElementList.hpp>
#include <Reaktoro/Core/SpeciesList.hpp>

namespace Reaktoro {

/// The class used to store and retrieve data of chemical species.
/// @see Element, Species
/// @ingroup Core
class Database
{
public:
    /// Return a Database object constructed with a given local file.
    /// @warning An exception is thrown if `path` does not point to a valid database file.
    /// @param path The path, including file name, to the database file.
    static auto fromFile(const String& path) -> Database;

    /// Return a Database object constructed with given database text contents.
    /// @param contents The contents of the database as a string.
    static auto fromContents(const String& contents) -> Database;

    /// Return a Database object constructed with given input stream containing the database text contents.
    /// @param stream The input stream containing the database file contents.
    static auto fromStream(std::istream& stream) -> Database;

    /// Construct a default Database object.
    Database();

    /// Construct a copy of a Database object.
    Database(const Database& other);

    /// Construct a Database object with given elements and species.
    Database(const Vec<Element>& elements, const Vec<Species>& species);

    /// Construct a Database object with given species (elements extracted from them).
    explicit Database(const Vec<Species>& species);

    /// Destroy this Database object.
    ~Database();

    /// Assign another Database object to this.
    auto operator=(Database other) -> Database&;

    /// Remove all species and elements from the database.
    auto clear() -> void;

    /// Add an element in the database.
    /// @note If an Element object with same symbol already exists in the
    /// Database container, the given Element object is not added.
    auto addElement(const Element& element) -> void;

    /// Add a species in the database.
    /// @note If a Species object with same name already exists in the Database
    /// container, the added Species object has its name slightly changed (e.g. `H2O` becomes `H2O!`).
    auto addSpecies(const Species& species) -> void;

    /// Add a list of species in the database.
    auto addSpecies(const Vec<Species>& species) -> void;

    /// Attach data to this database whose type is known at runtime only.
    auto attachData(const Any& data) -> void;

    /// Return all elements in the database.
    auto elements() const -> const ElementList&;

    /// Return all species in the database.
    auto species() const -> const SpeciesList&;

    /// Return all species in the database with given aggregate state.
    auto speciesWithAggregateState(AggregateState option) const -> SpeciesList;

    /// Return the attached data to this database whose type is known at runtime only.
    auto attachedData() const -> const Any&;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

} // namespace Reaktoro
