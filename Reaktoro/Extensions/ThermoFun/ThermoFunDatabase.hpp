// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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

/// The class used to store and retrieve data of chemical species from ThermoFun databases.
/// @ingroup ThermoFunExtension
class ThermoFunDatabase : public Database
{
public:
    /// Return a ThermoFunDatabase object initialized using an embedded database file.
    /// The currently supported database file names are:
    /// - `aq17`            (corresponding file: `aq17-thermofun.json`)
    /// - `cemdata18`       (corresponding file: `cemdata18-thermofun.json`)
    /// - `heracles`        (corresponding file: `heracles-thermofun.json`)
    /// - `mines16`         (corresponding file: `mines16-thermofun.json`)
    /// - `psinagra-12-07`  (corresponding file: `psinagra-12-07-thermofun.json`)
    /// - `slop98-organic`  (corresponding file: `slop98-organic-thermofun.json`)
    /// - `slop98`          (corresponding file: `slop98-thermofun.json`)
    /// @warning An exception is thrown if `name` is not one of the above names.
    /// @param name The name of the embedded ThermoFun database.
    static auto withName(const String& name) -> ThermoFunDatabase;

    /// Construct a ThermoFunDatabase object with given database file.
    /// If the string argument `file` does not have the contents of the
    /// database in JSON format or it is not a path to a valid local database
    /// file, an exception is thrown.
    /// @param file The contents or a path to a local file containing the database contents.
    static auto fromFile(const String& file) ->  ThermoFunDatabase;

    /// Construct a default ThermoFunDatabase object.
    ThermoFunDatabase();

    /// Construct a ThermoFunDatabase object with given name of embedded database file.
    /// For a list of currently supported names for embedded ThermoFun databases, see @ref ThermoFunDatabase::withName.
    explicit ThermoFunDatabase(const String& name);
};

} // namespace Reaktoro
