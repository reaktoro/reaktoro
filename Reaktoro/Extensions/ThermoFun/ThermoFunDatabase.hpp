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
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/Database.hpp>

// Forward declaration of ThermoFun classes
namespace ThermoFun { class Database; }

namespace Reaktoro {

/// The class used to store and retrieve data of chemical species from ThermoFun databases.
/// @ingroup ThermoFunExtension
class ThermoFunDatabase : public Database
{
public:
    /// Return a ThermoFunDatabase object constructed with an embedded database file.
    /// The currently supported database file names are:
    /// - `aq17`            (corresponding file: `aq17-thermofun.json`)
    /// - `cemdata18`       (corresponding file: `cemdata18-thermofun.json`)
    /// - `heracles`        (corresponding file: `heracles-thermofun.json`)
    /// - `mines16`         (corresponding file: `mines16-thermofun.json`)
    /// - `psinagra-12-07`  (corresponding file: `psinagra-12-07-thermofun.json`)
    /// - `slop98-organic`  (corresponding file: `slop98-organic-thermofun.json`)
    /// - `slop98`          (corresponding file: `slop98-thermofun.json`)
    /// @param name The name of the embedded ThermoFun database.
    /// @warning An exception is thrown if `name` is not one of the above names.
    static auto withName(const String& name) -> ThermoFunDatabase;

    /// Return a ThermoFunDatabase object constructed with a given local file.
    /// @param filepath The path, including file name, to the database file.
    /// @warning An exception is thrown if `filepath` does not point to a valid database file.
    static auto fromFile(const String& filepath) ->  ThermoFunDatabase;

    /// Return a ThermoFunDatabase object constructed with given local files.
    /// @param filepaths The paths, including file name, to the database files.
    /// @warning An exception is thrown if `filepaths` contain paths to invalid database files.
    static auto fromFiles(const Strings& filepaths) ->  ThermoFunDatabase;

    /// Return a ThermoFunDatabase object constructed with given database text contents.
    /// @param contents The contents of the database as a string.
    static auto fromContents(const String& contents) ->  ThermoFunDatabase;

    /// Disable all logs emitted by ThermoFun.
    static auto disableLogging() ->  void;

    /// Construct a default ThermoFunDatabase object.
    ThermoFunDatabase();

    /// Construct a ThermoFunDatabase object with given name of embedded database file.
    /// For a list of currently supported names for embedded ThermoFun databases, see @ref ThermoFunDatabase::withName.
    explicit ThermoFunDatabase(const String& name);

    /// Construct a ThermoFunDatabase object with given object of class ThermoFun::Database.
    explicit ThermoFunDatabase(const ThermoFun::Database& db);
};

} // namespace Reaktoro
