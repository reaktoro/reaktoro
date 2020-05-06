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
#include <Reaktoro/Core/Database.hpp>

namespace Reaktoro {

/// The class used to store and retrieve data of chemical species from PHREEQC databases.
class PhreeqcDatabase : public Database
{
public:
    /// Construct a default PhreeqcDatabase object.
    PhreeqcDatabase();

    /// Construct a PhreeqcDatabase object with given name of embedded PHREEQC database.
    /// This constructor initializes the PhreeqcDatabase object with an
    /// embedded PHREEQC database file. The following embedded file names are
    /// currently supported:"
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
    /// @param database The name of the embedded PHREEQC database file
    explicit PhreeqcDatabase(String database);

    /// Extend this PhreeqcDatabase object with contents in given database file.
    /// This method supports either a path to a database file, including its
    /// file name, or a multi-line string containing the database contents.
    /// @param database The path to the database file or its contents as a string
    auto load(String filename) -> PhreeqcDatabase&;

    /// Return the contents of an embedded PHREEQC database as a string.
    static auto contents(String database) -> String;
};

} // namespace Reaktoro
