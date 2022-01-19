// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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
#include <Reaktoro/Core/Database.hpp>

namespace Reaktoro {

/// The class used to store and retrieve data of chemical species from SUPCRT databases.
/// @ingroup Databases
class SupcrtDatabase : public Database
{
public:
    /// Construct a default SupcrtDatabase object.
    SupcrtDatabase();

    /// Construct a SupcrtDatabase object using an embedded database file.
    /// The currently supported database file names are:
    /// - `supcrt98`
    /// - `supcrt07`
    /// - `supcrt16`
    /// - `supcrtbl`
    /// @warning An exception is thrown if `name` is not one of the above names.
    /// @param name The name of the embedded database.
    SupcrtDatabase(const String& name);

    /// Construct a SupcrtDatabase object from another given Database object.
    /// @param other The Database object already containing chemical species data or empty.
    SupcrtDatabase(const Database& other);

    /// Return a SupcrtDatabase object initialized using an embedded database file.
    /// The currently supported database file names are:
    /// - `supcrt98`
    /// - `supcrt07`
    /// - `supcrt16`
    /// - `supcrtbl`
    /// @warning An exception is thrown if `name` is not one of the above names.
    /// @param name The name of the embedded database.
    static auto withName(const String& name) -> SupcrtDatabase;

    /// Return a SupcrtDatabase object constructed with a given local file.
    /// @warning An exception is thrown if `path` does not point to a valid database file.
    /// @param path The path, including file name, to the database file.
    static auto fromFile(const String& path) -> SupcrtDatabase;

    /// Return a SupcrtDatabase object initialized with given database text contents.
    /// @param contents The contents of the database as a string.
    static auto fromContents(const String& contents) -> SupcrtDatabase;
};

} // namespace Reaktoro
