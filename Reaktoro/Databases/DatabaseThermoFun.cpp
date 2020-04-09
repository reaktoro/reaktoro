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

#include "DatabaseThermoFun.hpp"

// Reaktoro includes
#include <Reaktoro/Core/Element.hpp>
#include <Reaktoro/Core/Species.hpp>

// ThermoFun includes
#include <ThermoFun/ThermoFun.h>

namespace Reaktoro {
namespace {

/// Convert a ThermoFun::Element object into a Reaktoro::Element object
auto convert(const ThermoFun::Element& element) -> Element
{
    Element converted;
    converted = converted.withName(element.name());
    converted = converted.withSymbol(element.symbol());
    converted = converted.withMolarMass(element.molarMass());
    return converted;
}

/// Return the elements and their coefficients in a species a ThermoFun::Substance object
/// @param db The ThermoFun::Database object
/// @param subtance The ThermoFun::Substance object
auto createElements(ThermoFun::Database& db, const ThermoFun::Substance& substance) -> Species::Elements
{
    Species::Elements elements;
    for(auto&& [element, coeff] : db.parseSubstanceFormula(substance.formula()))
    {
        if(element.symbol() == "Zz")
            continue;
        elements[convert(element)] = coeff;
    }
    return elements;
}

/// Convert a ThermoFun::Substance object into a Reaktoro::Species object
/// @param subtance The ThermoFun::Substance object
/// @param db The ThermoFun::Database object
auto convert(const ThermoFun::Substance& substance, ThermoFun::Database& db) -> Species
{
    Species species;
    species = species.withName(substance.symbol());
    species = species.withFormula(substance.formula());
    species = species.withElements(createElements(db, substance));
    species = species.withCharge(substance.charge());
    species = species.withAttachedData(substance);
    return species;
}

} // namespace

DatabaseThermoFun::DatabaseThermoFun()
{}

auto DatabaseThermoFun::fromFile(std::string path) ->  DatabaseThermoFun
{
    ThermoFun::Database thermofun_db(path);

    DatabaseThermoFun db;
    db.attachData(thermofun_db);

    for(auto [_, subs] : thermofun_db.mapSubstances())
        db.addSpecies(convert(subs, thermofun_db));

    return db;
}

} // namespace Reaktoro
