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

#include "EmbeddedDatabases.hpp"

// C++ includes
#include <map>

// Miniz includes
#include <miniz/zip_file.hpp>

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include "embedded/supcrt98.hpp"
#include "embedded/supcrt07.hpp"
#include "embedded/supcrt98-organics.hpp"
#include "embedded/supcrt07-organics.hpp"

namespace Reaktoro {
namespace {

/// The names of the built-in SUPCRT databases
const std::vector<std::string> supcrt_databases =
{
    "supcrt98.xml",
    "supcrt07.xml",
    "supcrt98-organics.xml",
    "supcrt07-organics.xml",
};

/// The arrays containing the zipped data of the built-in SUPCRT databases
const std::vector<unsigned char*> supcrt_databases_data =
{
    supcrt98_zip,
    supcrt07_zip,
    supcrt98_organics_zip,
    supcrt07_organics_zip,
};

/// The lengths of the arrays containing the zipped data of the built-in SUPCRT databases
const std::vector<unsigned int> supcrt_databases_len =
{
    supcrt98_zip_len,
    supcrt07_zip_len,
    supcrt98_organics_zip_len,
    supcrt07_organics_zip_len,
};

} // namespace

auto supcrtEmbeddedDatabaseTextContent(std::string name) -> std::string
{
    // Get the index of the database, either named, e.g., supcrt98.xml or supcrt98
    const auto i = index(supcrt_databases, name);
    const auto j = index(supcrt_databases, name + ".xml");
    const auto idx = std::min(i, j);

    // Return empty string if there is no built-in database if such name
    if(idx >= supcrt_databases.size())
        return "";

    // If the database name was provided without extension, then add it
    if(j < supcrt_databases.size())
        name += ".xml";

    // The begin and end pointers to the array containing the database data
    const auto& begin = supcrt_databases_data[idx];
    const auto& end = begin + supcrt_databases_len[idx];

    // Create the vector of unsigned char
    std::vector<unsigned char> data(begin, end);

    // Read the zip file
    zip_file file(data);

    // Return the contents of the unzipped file as a string
    return file.read(name);
}

auto supcrtEmbeddedDatabaseNames() -> std::vector<std::string>
{
    return supcrt_databases;
}

} // namespace Reaktoro
