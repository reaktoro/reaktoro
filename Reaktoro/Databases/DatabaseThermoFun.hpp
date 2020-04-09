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

// C++ includes
#include <string>

// Reaktoro includes
#include <Reaktoro/Core/Database.hpp>

namespace Reaktoro {

/// The class used to store and retrieve data of chemical species from ThermoFun databases.
/// @ingroup Databases
class DatabaseThermoFun : public Database
{
public:
    /// Construct a default DatabaseThermoFun object.
    DatabaseThermoFun();

    /// Construct a DatabaseThermoFun object with given database file.
    /// If `path` does not point to a valid database file, an exception is thrown.
    /// @param path The path, including file name, to the database file.
    static auto fromFile(std::string path) ->  DatabaseThermoFun;
};

} // namespace Reaktoro
