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

#include "NasaDatabase.hpp"

// C++ includes
#include <fstream>

// CMakeRC includes
#include <cmrc/cmrc.hpp>

CMRC_DECLARE(ReaktoroDatabases);

// Reaktoro includes
#include <Reaktoro/Extensions/Nasa/NasaDatabaseParseUtils.hpp>

namespace Reaktoro {

/// Return the contents of the embedded NASA database with given name (or empty)
auto getNasaDatabaseContent(String name) -> String
{
    error(!oneof(name,
        "thermo.inp"),
        "Could not load embedded NASA database file with name `", name, "`. ",
        "The currently supported names are: \n"
        "    - thermo.inp    \n",
        "");
    auto fs = cmrc::ReaktoroDatabases::get_filesystem();
    auto contents = fs.open("databases/nasa/" + name);
    return String(contents.begin(), contents.end());
}

NasaDatabase::NasaDatabase()
: Database()
{}

NasaDatabase::NasaDatabase(String name)
: NasaDatabase(NasaDatabase::withName(name))
{}

auto NasaDatabase::withName(String name) -> NasaDatabase
{
    const auto content = getNasaDatabaseContent(name);

    NasaDatabase db;
    // detail::NasaDatabaseHelper helper(content);
    // db.addSpecies(helper.species);
    // db.attachData(helper);
    return db;
}

auto NasaDatabase::fromFile(String path) -> NasaDatabase
{
    std::ifstream file(path);

    error(!file.is_open(), "Could not open NASA database file with given path: ", path);

    const auto lines = NasaUtils::createTextLines(file);
    const auto lines_prods = NasaUtils::getTextLinesForProducts(lines);
    const auto lines_reacs = NasaUtils::getTextLinesForReactants(lines);

    const auto species_prods = NasaUtils::createNasaSpeciesVector(lines_prods);
    const auto species_reacs = NasaUtils::createNasaSpeciesVector(lines_reacs);

    NasaDatabase db;
    // db.addSpecies(convertNasaSpeciesVector(species_prods));
    // db.addSpecies(convertNasaSpeciesVector(species_reacs));

    file.close();

    return db;
}

} // namespace Reaktoro
