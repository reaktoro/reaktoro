// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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

#include "ThermoFunDatabase.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/Element.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Extensions/ThermoFun/ThermoFunEngine.hpp>

// ThermoFun includes
#include <ThermoFun/ThermoFun.h>

namespace Reaktoro {
namespace {

/// Return the standard thermodynamic property function of a species with given name.
auto createStandardThermoModel(const ThermoFunEngine& engine, const String& species) -> StandardThermoModel
{
    return [=](real T, real P) -> StandardThermoProps
    {
        return engine.props(T, P, species);
    };
}

/// Convert a ThermoFun::Element object into a Reaktoro::Element object
auto createElement(const ThermoFun::Element& element) -> Element
{
    Element converted;
    converted = converted.withName(element.name());
    converted = converted.withSymbol(element.symbol());
    converted = converted.withMolarMass(element.molarMass() * 1e-3); // from g/mol to kg/mol
    return converted;
}

/// Return the elements and their coefficients in a species a ThermoFun::Substance object
auto createElements(const ThermoFunEngine& engine, const ThermoFun::Substance& substance) -> Map<Element, double>
{
    auto db = engine.database();
    Map<Element, double> elements;
    for(auto&& [element, coeff] : db.parseSubstanceFormula(substance.formula()))
    {
        if(element.symbol() == "Zz")
            continue; // skip charge element in ThermoFun (charge is explicitly accessed in Reaktoro::Species)
        elements[createElement(element)] = coeff;
    }
    return elements;
}

/// Convert a ThermoFun::Substance object into a Reaktoro::Species object
auto createSpecies(const ThermoFunEngine& engine, const ThermoFun::Substance& substance) -> Species
{
    Species species;
    species = species.withName(substance.symbol());
    species = species.withFormula(substance.formula());
    species = species.withSubstance(substance.name());
    species = species.withElements(createElements(engine, substance));
    species = species.withCharge(substance.charge());
    // species = species.withAggregateState(aggregateState(substance));
    species = species.withStandardThermoModel(createStandardThermoModel(engine, substance.symbol()));
    species = species.withAttachedData(substance);
    return species;
}

} // namespace

ThermoFunDatabase::ThermoFunDatabase()
{}

auto ThermoFunDatabase::fromFile(const String& path) ->  ThermoFunDatabase
{
    ThermoFun::Database database(path);

    ThermoFunEngine engine(database);

    ThermoFunDatabase db;
    db.attachData(engine);
    for(auto [_, subs] : database.mapSubstances())
        db.addSpecies(createSpecies(engine, subs));

    return db;
}

} // namespace Reaktoro
