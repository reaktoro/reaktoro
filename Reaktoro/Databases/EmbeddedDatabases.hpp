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
#include <string>
#include <vector>

namespace Reaktoro {

/// Return the text content of an embedded SUPCRT database as a string.
/// This method searches among all built-in SUPCRT databases and
/// return the text content of the database file as a string.
/// If the given database `name` is not found, an empty string is returned.
/// @param name The name of the embedded SUPCRT database.
auto supcrtEmbeddedDatabaseTextContent(std::string name) -> std::string;

/// Return the list of names of all embedded SUPCRT databases.
auto supcrtEmbeddedDatabaseNames() -> std::vector<std::string>;

} // namespace Reaktoro
