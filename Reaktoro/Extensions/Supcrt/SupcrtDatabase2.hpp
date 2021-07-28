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
class SupcrtDatabase2 : public Database
{
public:
    /// Construct a default SupcrtDatabase2 object.
    SupcrtDatabase2();

    /// Construct a SupcrtDatabase2 object using an embedded database file.
    /// The currently supported database file names are:
    /// - `supcrt98`
    /// - `supcrt07`
    /// - `supcrt16`
    /// - `supcrtbl`
    /// @warning An exception is thrown if `name` is not one of the above names.
    /// @param name The name of the embedded database.
    SupcrtDatabase2(const String& name);

    /// Construct a SupcrtDatabase2 object from another given Database object.
    /// @param other The Database object already containing chemical species data or empty.
    SupcrtDatabase2(const Database& other);

    /// Return a SupcrtDatabase2 object initialized using an embedded database file.
    /// The currently supported database file names are:
    /// - `supcrt98`
    /// - `supcrt07`
    /// - `supcrt16`
    /// - `supcrtbl`
    /// @warning An exception is thrown if `name` is not one of the above names.
    /// @param name The name of the embedded database.
    static auto withName(const String& name) -> SupcrtDatabase2;

    /// Return a SupcrtDatabase2 object initialized using a given path to a database file.
    /// @warning An exception is thrown if `path` does not point to a valid database file.
    /// @param path The path, including file name, to the database file.
    static auto fromFile(const String& path) -> SupcrtDatabase2;
};

} // namespace Reaktoro
