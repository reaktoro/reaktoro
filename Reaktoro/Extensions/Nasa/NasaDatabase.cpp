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

#include "NasaDatabase.hpp"

// Reaktoro includes
#include <Reaktoro/Core/Embedded.hpp>

namespace Reaktoro {

NasaDatabase::NasaDatabase(Database database)
: Database(std::move(database))
{}

NasaDatabase::NasaDatabase(String name)
: NasaDatabase(NasaDatabase::withName(name))
{}

auto NasaDatabase::withName(String name) -> NasaDatabase
{
    errorif(!oneof(name, "nasa-cea"),
        "Could not load embedded NASA database file with name `", name, "`. ",
        "The currently supported names are: \n"
        "    - nasa-cea \n",
        "");
    const String contents = Embedded::get("databases/reaktoro/" + name + ".json");
    return fromStringJSON(contents);
}

} // namespace Reaktoro
