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
    /// Construct a default Database instance.
    Database();

    /// Construct a copy of a Database instance.
    Database(const Database& other);

    /// Construct a Database instance with given species.
    explicit Database(const Vec<Species>& species);

    /// Destroy this Database instance.
    ~Database();

    /// Assign another Database instance to this.
    auto operator=(Database other) -> Database&;

    /// Remove all species and elements from the database.
    auto clear() -> void;

    /// Add a species in the database.
    auto addSpecies(const Species& species) -> void;

    /// Add a list of species in the database.
    auto addSpecies(const Vec<Species>& species) -> void;

    /// Attach data to this database whose type is known at runtime only.
    auto attachData(const Any& data) -> void;

    /// Return all elements in the database.
    auto elements() const -> ElementListConstRef;

    /// Return all species in the database.
    auto species() const -> SpeciesListConstRef;

    /// Return all species in the database with given aggregate state.
    auto speciesWithAggregateState(AggregateState option) const -> SpeciesList;

    /// Return the attached data to this database whose type is known at runtime only.
    auto attachedData() const -> const Any&;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
