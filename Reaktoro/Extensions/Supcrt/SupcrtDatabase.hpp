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
    /// If `name` does not correspond to one of the following names, an exception is thrown:
    /// - `supcrt98`
    /// - `supcrt07`
    /// - `supcrt98-organics`
    /// - `supcrt07-organics`
    /// @param name The name of the embedded SUPCRT database.
    SupcrtDatabase(String name);

    /// Return a SupcrtDatabase object initialized using an embedded database file.
    /// If `name` does not correspond to one of the following names, an exception is thrown:
    /// - `supcrt98`
    /// - `supcrt07`
    /// - `supcrt98-organics`
    /// - `supcrt07-organics`
    /// @param name The name of the embedded SUPCRT database.
    static auto withName(String name) -> SupcrtDatabase;

    /// Return a SupcrtDatabase object initialized using a given path to a database file.
    /// If `path` does not point to a valid database file, an exception is thrown.
    /// @param path The path, including file name, to the database file.
    static auto fromFile(String path) -> SupcrtDatabase;
};

} // namespace Reaktoro
