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

#include "AggregateState.hpp"

// Reaktoro includes
#include <Reaktoro/Core/ChemicalFormula.hpp>

namespace Reaktoro {

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

auto identifyAggregateState(const std::string& name) -> AggregateState
{
    auto contains = [&](auto str) { return name.rfind(str) != std::string::npos; };
    auto charged = [&](auto str) { return parseElectricCharge(name) != 0; };
    if(contains("(g)"))   return AggregateState::Gas;
    if(contains("(l)"))   return AggregateState::Liquid;
    if(contains("(s)"))   return AggregateState::Solid;
    if(contains("(pl)"))  return AggregateState::Plasma;
    if(contains("(cd)"))  return AggregateState::CondensedPhase;
    if(contains("(fl)"))  return AggregateState::Fluid;
    if(contains("(lc)"))  return AggregateState::LiquidCrystal;
    if(contains("(cr)"))  return AggregateState::CrystallineSolid;
    if(contains("(am)"))  return AggregateState::AmorphousSolid;
    if(contains("(vit)")) return AggregateState::Vitreous;
    if(contains("(ads)")) return AggregateState::Adsorbed;
    if(contains("(mon)")) return AggregateState::Monomeric;
    if(contains("(pol)")) return AggregateState::Polymeric;
    if(contains("(ss)"))  return AggregateState::SolidSolution;
    if(contains("(ex)"))  return AggregateState::IonExchange;
    if(contains("(aq)"))  return AggregateState::Aqueous;
    if(charged(name)) return AggregateState::Aqueous;
    return AggregateState::Undefined;
}

} // namespace Reaktoro
