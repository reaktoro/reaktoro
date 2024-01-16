// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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
#include <Reaktoro/Core/Database.hpp>

// Reaktoro includes
namespace Reaktoro {

/// Used to support thermodynamic databases in NASA CEA format.
/// @ingroup Databases
class NasaDatabase : public Database
{
public:
    // Inherit all Database constructors.
    using Database::Database;

    /// Construct a NasaDatabase object using an object of Database.
    NasaDatabase(Database database);

    /// Construct a NasaDatabase object using an embedded database file.
    /// The currently supported embedded NASA database files are named:
    /// - `nasa-cea`
    /// @param name The name of the embedded NASA database file
    /// @warning An exception is thrown if `name` is not one of the above supported names.
    NasaDatabase(String name);

    /// Return a NasaDatabase object initialized using an embedded NASA database.
    /// The currently supported embedded NASA database files are named:
    /// - `nasa-cea`
    /// @param name The name of the embedded NASA database file
    /// @warning An exception is thrown if `name` is not one of the above supported names.
    static auto withName(String name) -> NasaDatabase;
};

} // namespace Reaktoro
