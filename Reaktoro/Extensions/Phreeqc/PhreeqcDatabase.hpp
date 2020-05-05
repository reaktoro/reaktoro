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

    /// Construct a PhreeqcDatabase object with given database.
    /// This constructor supports initialization of the PhreeqcDatabase object
    /// with either a path to the database file, including its file name, or a
    /// multi-line string containing the database contents itself.
    /// @param database The path to the database file or its contents in a string
    explicit PhreeqcDatabase(String database);

    /// Extend this PhreeqcDatabase object with contents in given database file.
    /// This method supports either a path to a database file, including its
    /// file name, or a multi-line string containing the database contents.
    /// @param database The path to the database file or its contents as a string
    auto load(String filename) -> PhreeqcDatabase&;
};

} // namespace Reaktoro
