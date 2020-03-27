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

#include "DatabaseThermoFun.hpp"

// C++ includes

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/Element.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Databases/EmbeddedDatabases.hpp>

// ThermoFun includes
#include <ThermoFun/ThermoFun.h>

namespace Reaktoro {
namespace {

/// Convert a ThermoFun::Element symbol in a Reaktoro::Element name.
auto convertThermoFunElementSymbol(const ThermoFun::Element& element) -> std::string
{
    /// ThermoFun::Element::symbol is equivalent to Reaktoro::Element::name except for electric charge, which is Z instead of Zz
    return element.symbol() == "Zz" ? "Z" : element.symbol();
}

/// Return the elements and their coefficients in a species a ThermoFun::Substance object
/// @param db The ThermoFun::Database object
/// @param subtance The ThermoFun::Substance object
/// @param element_map The unique and immutable Reaktoro::Element objects previously created
auto speciesElements(ThermoFun::Database& db, const ThermoFun::Substance& substance, const ElementMap& element_map) -> std::map<Element, double>
{
    std::map<Element, double> elements;

    for(auto [symbol, coeff] : db.parseSubstanceFormula(substance.formula()))
    {
        const auto name = convertThermoFunElementSymbol(symbol);
        const auto element = element_map.at(name);
        elements[element] = coeff;
    }

    return elements;
}

/// Convert a ThermoFun::Element object into a Reaktoro::Element object
auto convertThermoFunElement(const ThermoFun::Element& element) -> Element
{
    Element converted;
    converted = converted.withName(convertThermoFunElementSymbol(element));
    converted = converted.withMolarMass(element.molarMass());
    return converted;
}

/// Convert a ThermoFun::Substance object into a Reaktoro::Species object
/// @param subtance The ThermoFun::Substance object
/// @param db The ThermoFun::Database object
/// @param element_map The unique and immutable Reaktoro::Element objects previously created
auto convertThermoFunSubstance(const ThermoFun::Substance& substance, ThermoFun::Database& db, const ElementMap& element_map) -> Species
{
    Species species;
    species = species.withName(substance.symbol());
    species = species.withFormula(substance.formula());
    species = species.withElements(speciesElements(db, substance, element_map));
    species = species.withData(substance);
    return species;
}

/// Create the Reaktoro::Element objects for the Reaktoro::Database object
auto createElementMap(ThermoFun::Database& db) -> ElementMap
{
    ElementMap element_map;

    for(auto [_, elem] : db.mapElements())
    {
        Element element = convertThermoFunElement(elem);
        element_map[element.name()] = element;
    }

    return element_map;
}

/// Create the Reaktoro::Species objects for the Reaktoro::Database object
auto createSpeciesMap(ThermoFun::Database& db, const ElementMap& element_map) -> SpeciesMap
{
    SpeciesMap species_map;

    for(auto [_, subs] : db.mapSubstances())
    {
        Species species = convertThermoFunSubstance(subs, db, element_map);
        species_map[species.name()] = species;
    }

    return species_map;
}

} // namespace

DatabaseThermoFun::DatabaseThermoFun()
{}

auto DatabaseThermoFun::fromFile(std::string path) ->  DatabaseThermoFun
{
    ThermoFun::Database thermofun_db(path);

    ElementMap element_map = createElementMap(thermofun_db);
    SpeciesMap species_map = createSpeciesMap(thermofun_db, element_map);

    DatabaseThermoFun db;
    db.setElements(element_map);
    db.setSpecies(species_map);
    db.setData(thermofun_db);
    return db;
}

} // namespace Reaktoro
