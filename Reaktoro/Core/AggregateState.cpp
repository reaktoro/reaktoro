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

#include "AggregateState.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/NamingUtils.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Core/ChemicalFormula.hpp>

namespace Reaktoro {

auto operator<<(std::ostream& out, AggregateState option) -> std::ostream&
{
    switch(option)
    {
    case AggregateState::Gas:              out << "Gas";              return out;
    case AggregateState::Liquid:           out << "Liquid";           return out;
    case AggregateState::Solid:            out << "Solid";            return out;
    case AggregateState::Plasma:           out << "Plasma";           return out;
    case AggregateState::CondensedPhase:   out << "CondensedPhase";   return out;
    case AggregateState::Fluid:            out << "Fluid";            return out;
    case AggregateState::LiquidCrystal:    out << "LiquidCrystal";    return out;
    case AggregateState::CrystallineSolid: out << "CrystallineSolid"; return out;
    case AggregateState::AmorphousSolid:   out << "AmorphousSolid";   return out;
    case AggregateState::Vitreous:         out << "Vitreous";         return out;
    case AggregateState::Adsorbed:         out << "Adsorbed";         return out;
    case AggregateState::Monomeric:        out << "Monomeric";        return out;
    case AggregateState::Polymeric:        out << "Polymeric";        return out;
    case AggregateState::SolidSolution:    out << "SolidSolution";    return out;
    case AggregateState::IonExchange:      out << "IonExchange";      return out;
    case AggregateState::Aqueous:          out << "Aqueous";          return out;
    default:                               out << "Undefined";        return out;
    }
}

auto parseAggregateState(const std::string& symbol) -> AggregateState
{
    if(symbol == "g")   return AggregateState::Gas;
    if(symbol == "l")   return AggregateState::Liquid;
    if(symbol == "s")   return AggregateState::Solid;
    if(symbol == "pl")  return AggregateState::Plasma;
    if(symbol == "cd")  return AggregateState::CondensedPhase;
    if(symbol == "fl")  return AggregateState::Fluid;
    if(symbol == "lc")  return AggregateState::LiquidCrystal;
    if(symbol == "cr")  return AggregateState::CrystallineSolid;
    if(symbol == "am")  return AggregateState::AmorphousSolid;
    if(symbol == "vit") return AggregateState::Vitreous;
    if(symbol == "ads") return AggregateState::Adsorbed;
    if(symbol == "mon") return AggregateState::Monomeric;
    if(symbol == "pol") return AggregateState::Polymeric;
    if(symbol == "ss")  return AggregateState::SolidSolution;
    if(symbol == "ex")  return AggregateState::IonExchange;
    if(symbol == "aq")  return AggregateState::Aqueous;
    return AggregateState::Undefined;
}

auto identifyAggregateState(const std::string& substance) -> AggregateState
{
    const auto [name, suffix] = splitSpeciesNameSuffix(substance);
    const auto words = split(suffix, " ,"); // suffix = "s, calcite" => words = {"s", "calcite"}
    auto charged = [&](auto str) { return parseElectricCharge(str) != 0; };
    if(contains(words, "g"))   return AggregateState::Gas;
    if(contains(words, "l"))   return AggregateState::Liquid;
    if(contains(words, "s"))   return AggregateState::Solid;
    if(contains(words, "pl"))  return AggregateState::Plasma;
    if(contains(words, "cd"))  return AggregateState::CondensedPhase;
    if(contains(words, "fl"))  return AggregateState::Fluid;
    if(contains(words, "lc"))  return AggregateState::LiquidCrystal;
    if(contains(words, "cr"))  return AggregateState::CrystallineSolid;
    if(contains(words, "am"))  return AggregateState::AmorphousSolid;
    if(contains(words, "vit")) return AggregateState::Vitreous;
    if(contains(words, "ads")) return AggregateState::Adsorbed;
    if(contains(words, "mon")) return AggregateState::Monomeric;
    if(contains(words, "pol")) return AggregateState::Polymeric;
    if(contains(words, "ss"))  return AggregateState::SolidSolution;
    if(contains(words, "ex"))  return AggregateState::IonExchange;
    if(contains(words, "aq"))  return AggregateState::Aqueous;
    if(charged(name)) return AggregateState::Aqueous;
    return AggregateState::Undefined;
}

} // namespace Reaktoro
