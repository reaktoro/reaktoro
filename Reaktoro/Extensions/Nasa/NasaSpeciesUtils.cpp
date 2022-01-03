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

#include "NasaSpeciesUtils.hpp"

// Reaktoro includes
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Extensions/Nasa/NasaSpecies.hpp>
#include <Reaktoro/Extensions/Nasa/NasaThermoModels.hpp>
#include <Reaktoro/Singletons/Elements.hpp>

namespace Reaktoro {
namespace NasaUtils {

auto correctElementSymbol(const String& symbol) -> String
{
    String corrected = lowercase(symbol);
    corrected[0] = std::toupper(corrected[0]);
    return corrected;
}

auto createElement(const String& symbol0) -> Element
{
    const String symbol = correctElementSymbol(symbol0);
    const auto element = Elements::withSymbol(symbol);
    error(!element.has_value(), "Could not create Element object with given NASA CEA element symbol: ", symbol0);
    return element.value();
}

auto createElements(const NasaSpecies& species) -> ElementalComposition
{
    Pairs<Element, double> pairs;
    for(auto const& [symbol, coeff] : species.formula)
        if(symbol != "E") // skip charge, which is handled with Species::withCharge and Species::charge methods
            pairs.emplace_back(createElement(symbol), coeff);
    return ElementalComposition(pairs);
}

auto charge(const NasaSpecies& species) -> double
{
    for(auto const& [symbol, coeff] : species.formula)
        if(symbol == "E")
            return -coeff;
    return 0.0;
}

auto aggregateState(const NasaSpecies& species) -> AggregateState
{
    return (species.aggregatestate == NasaAggregateState::Gas) ?
        AggregateState::Gas :
        AggregateState::CondensedPhase;
}

auto tags(const NasaSpecies& species) -> Strings
{
    if(species.type == NasaSpeciesType::Product)
        return {"product"};
    if(species.type == NasaSpeciesType::Reactant)
        return {"reactant"};
    if(species.type == NasaSpeciesType::ProductReactant)
        return {"product", "reactant"};
    return {};
}

auto convertSpecies(const NasaSpecies& species) -> Species
{
    return Species()
        .withName(species.name)
        .withFormula(species.name)
        .withSubstance(species.name)
        .withElements(NasaUtils::createElements(species))
        .withCharge(NasaUtils::charge(species))
        .withAggregateState(NasaUtils::aggregateState(species))
        .withStandardThermoModel(StandardThermoModelNasa(species))
        .withTags(NasaUtils::tags(species))
        .withAttachedData(species)
        ;
}

} // namespace NasaUtils
} // namespace Reaktoro
