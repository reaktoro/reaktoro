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
#include <map>
#include <memory>
#include <string>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Core/Database.hpp>

namespace Reaktoro {

/// The class used to store and retrieve data of chemical species from SUPCRT databases.
/// @see Element, Species
/// @ingroup Databases
class DatabaseSupcrt : public Database
{
public:
    /// Construct a default DatabaseSupcrt object.
    DatabaseSupcrt();

    /// Construct a DatabaseSupcrt object with given name of a built-in database file.
    /// If `name` does not correspond to one of the following names, an exception is thrown:
    /// - `supcrt98`
    /// - `supcrt07`
    /// - `supcrt98-organics`
    /// - `supcrt07-organics`
    /// @param name The name of the built-in SUPCRT database.
    static auto withName(std::string name) ->  DatabaseSupcrt;

    /// Construct a DatabaseSupcrt object with given database file.
    /// If `path` does not point to a valid database file, an exception is thrown.
    /// @param path The path, including file name, to the database file.
    static auto fromFile(std::string path) ->  DatabaseSupcrt;
};

} // namespace Reaktoro
