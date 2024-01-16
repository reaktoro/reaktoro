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

// Forward declaration of PHREEQC class
class PHREEQC;

namespace Reaktoro {

/// The class used to store and retrieve data of chemical species from PHREEQC databases.
class PhreeqcDatabase : public Database
{
public:
    /// Construct a default PhreeqcDatabase object.
    PhreeqcDatabase();

    /// Construct a PhreeqcDatabase object using an embedded PHREEQC database.
    /// If `name` does not correspond to one of the following names, an exception is thrown:
    /// - Amm.dat
    /// - frezchem.dat
    /// - iso.dat
    /// - llnl.dat
    /// - minteq.dat
    /// - minteq.v4.dat
    /// - phreeqc.dat
    /// - pitzer.dat
    /// - sit.dat
    /// - wateq4f.dat
    /// @param name The name of the embedded PHREEQC database file
    PhreeqcDatabase(const String& name);

    /// Extend this PhreeqcDatabase object with contents in given database file.
    /// This method supports either a path to a database file, including its
    /// file name, or a multi-line string containing the database contents.
    /// @param database The path to the database file or its contents as a string
    auto load(const String& database) -> PhreeqcDatabase&;

    /// Return the underlying PHREEQC object.
    auto ptr() const -> SharedPtr<PHREEQC> const&;

    /// Return a PhreeqcDatabase object constructed with an embedded database file.
    /// If `name` does not correspond to one of the following names, an exception is thrown:
    /// - Amm.dat
    /// - frezchem.dat
    /// - iso.dat
    /// - llnl.dat
    /// - minteq.dat
    /// - minteq.v4.dat
    /// - phreeqc.dat
    /// - pitzer.dat
    /// - sit.dat
    /// - wateq4f.dat
    /// @param name The name of the embedded PHREEQC database file
    static auto withName(const String& name) -> PhreeqcDatabase;

    /// Return a PhreeqcDatabase object constructed with a given local file.
    /// @param path The path, including file name, to the database file.
    /// @warning An exception is thrown if `path` does not point to a valid database file.
    static auto fromFile(const String& path) -> PhreeqcDatabase;

    /// Return a PhreeqcDatabase object constructed with given database text contents.
    /// @param contents The contents of the database as a string.
    static auto fromContents(const String& contents) -> PhreeqcDatabase;

    /// Return the contents of an embedded PHREEQC database as a string.
    /// @param name The name of the embedded PHREEQC database file
    /// @see withName
    static auto contents(const String& name) -> String;

    /// Return the names of the currently supported embedded PHREEQC databases.
    static auto namesEmbeddedDatabases() -> Strings;

private:
    /// The underlying PHREEQC object containing the state of PHREEQC after parsing the database.
    SharedPtr<PHREEQC> m_ptr;
};

} // namespace Reaktoro
