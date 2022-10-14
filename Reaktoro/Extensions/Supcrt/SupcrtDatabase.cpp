// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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

#include "SupcrtDatabase.hpp"

// C++ includes
#include <fstream>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/Data.hpp>
#include <Reaktoro/Core/Embedded.hpp>
#include <Reaktoro/Core/Support/DatabaseParser.hpp>

namespace Reaktoro {

SupcrtDatabase::SupcrtDatabase()
: Database()
{}

SupcrtDatabase::SupcrtDatabase(const Database& other)
: Database(other)
{}

SupcrtDatabase::SupcrtDatabase(const String& name)
: Database(SupcrtDatabase::withName(name))
{}

auto SupcrtDatabase::withName(const String& name) -> SupcrtDatabase
{
    errorif(!oneof(name,
        "supcrt98",
        "supcrt07",
        "supcrt16",
        "supcrtbl",
        "supcrtbl-organics"),
        "Could not load embedded database file with name `", name, "`. ",
        "The currently supported names are: \n"
        "    - supcrt98 \n",
        "    - supcrt07 \n",
        "    - supcrt16 \n",
        "    - supcrtbl \n",
        "    - supcrtbl-organics \n",
        "");
    const auto text = Embedded::get("databases/reaktoro/" + name + ".yaml");
    const auto doc = Data::parse(text);
    DatabaseParser dbparser(doc);
    return Database(dbparser);
}

auto SupcrtDatabase::fromFile(const String& path) -> SupcrtDatabase
{
    std::ifstream file(path);
    errorif(!file.is_open(),
        "Could not open file `", path, "`. Ensure the given file path "
        "is relative to the directory where your application is RUNNING "
        "(not necessarily where the executable is located!). Alternatively, "
        "try a full path to the file (e.g., "
        "in Windows, `C:\\User\\username\\mydata\\mydatabase.yaml`, "
        "in Linux and macOS, `/home/username/mydata/mydatabase.yaml`).");
    auto doc = Data::parse(file);
    DatabaseParser dbparser(doc);
    return Database(dbparser);
}

auto SupcrtDatabase::fromContents(const String& contents) -> SupcrtDatabase
{
    auto doc = Data::parse(contents);
    DatabaseParser dbparser(doc);
    return Database(dbparser);
}

} // namespace Reaktoro
