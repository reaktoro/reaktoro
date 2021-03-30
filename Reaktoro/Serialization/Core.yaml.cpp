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

#include "Core.yaml.hpp"

// Reaktoro includes
#include <Reaktoro/Common/YAML.hpp>
#include <Reaktoro/Core/ChemicalFormula.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Element.hpp>
#include <Reaktoro/Core/FormationReaction.hpp>
#include <Reaktoro/Core/Param.hpp>
#include <Reaktoro/Core/Params.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Serialization/Common.yaml.hpp>

namespace Reaktoro {

//=====================================================================================================================

auto operator<<(yaml& node, const AggregateState& obj) -> void
{
    std::stringstream ss;
    ss << obj;
    node = ss.str();
}

auto operator>>(const yaml& node, AggregateState& obj) -> void
{
    obj = parseAggregateState(node.as<std::string>());
}

//=====================================================================================================================

auto operator<<(yaml& node, const ChemicalFormula& obj) -> void
{
    node = obj.str();
}

auto operator>>(const yaml& node, ChemicalFormula& obj) -> void
{
    obj = ChemicalFormula(node.as<std::string>());
}

//=====================================================================================================================

auto operator<<(yaml& node, const ChemicalSystem& obj) -> void
{
    errorif(true, "`auto operator<<(yaml& node, const ChemicalSystem& obj) -> void` not implemented!");
}

auto operator>>(const yaml& node, ChemicalSystem& obj) -> void
{
    errorif(true,  "`auto operator<<=(ChemicalSystem& obj) -> void` not implemented!");
}

//=====================================================================================================================

auto operator<<(yaml& node, const Element& obj) -> void
{
    node["Symbol"]            = obj.symbol();
    node["Name"]              = obj.name();
    node["AtomicNumber"]      = obj.atomicNumber();
    node["AtomicWeight"]      = obj.atomicWeight();
    node["Electronegativity"] = obj.electronegativity();
    node["Tags"]              = obj.tags();
}

auto operator>>(const yaml& node, Element& obj) -> void
{
    Element::Attribs attribs;
    node.required("Symbol", attribs.symbol);
    node.required("Name", attribs.name);
    node.required("AtomicNumber", attribs.atomic_number);
    node.required("AtomicWeight", attribs.atomic_weight);
    node.required("Electronegativity", attribs.electronegativity);
    node.required("Tags", attribs.tags);
    obj = Element(attribs);
}

//=====================================================================================================================

auto operator<<(yaml& node, const ElementalComposition& obj) -> void
{
    errorif(true, "`auto operator<<(yaml& node, const ElementalComposition& obj) -> void` not implemented!");
}

auto operator>>(const yaml& node, ElementalComposition& obj) -> void
{
    errorif(true, "`auto operator>>(const yaml& node, ElementalComposition& obj) -> void` not implemented!");
}

//=====================================================================================================================

auto operator<<(yaml& node, const FormationReaction& obj) -> void
{
    errorif(true, "`auto operator<<(yaml& node, const FormationReaction& obj) -> void` not implemented!");
}

auto operator>>(const yaml& node, FormationReaction& obj) -> void
{
    errorif(true, "`auto operator>>(const yaml& node, FormationReaction& obj) -> void` not implemented!");
}

//=====================================================================================================================

auto operator<<(yaml& node, const Param& obj) -> void
{
    node = obj.value();
}

auto operator>>(const yaml& node, Param& obj) -> void
{
    obj = node.as<double>();
}

//=====================================================================================================================

auto operator<<(yaml& node, const Params& obj) -> void
{
    node = obj.data();
}

auto operator>>(const yaml& node, Params& obj) -> void
{
    auto values = node.as<Vec<double>>();
    obj = Params(values.begin(), values.end());
}

//=====================================================================================================================

auto operator<<(yaml& node, const Phase& obj) -> void
{
    errorif(true, "`auto operator<<(yaml& node, const Phase& obj) -> void` not implemented!");
}

auto operator>>(const yaml& node, Phase& obj) -> void
{
    errorif(true,  "`auto operator<<=(Phase& obj) -> void` not implemented!");
}

//=====================================================================================================================

auto operator<<(yaml& node, const Species& obj) -> void
{
    errorif(true, "`auto operator<<(yaml& node, const Species& obj) -> void` not implemented!");
}

auto operator>>(const yaml& node, Species& obj) -> void
{
    Species::Attribs attribs;
    node.optional("Name", attribs.name);
    node.required("Formula", attribs.formula);
    node.optional("Substance", attribs.substance);
    node.optional("Elements", attribs.elements);
    node.optional("Charge", attribs.charge);
    node.optional("AggregateState", attribs.aggregate_state);
    node.optional("FormationReaction", attribs.formation_reaction);
    // node.optional("StandardThermoModel", attribs.std_thermo_model);
    node.optional("Tags", attribs.tags);
    obj = Species(attribs);
}

} // namespace YAML
