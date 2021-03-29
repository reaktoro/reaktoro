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
#include <Reaktoro/Core/Param.hpp>
#include <Reaktoro/Core/Params.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Serialization/Common.yaml.hpp>

namespace Reaktoro {

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
    node["MolarMass"]         = obj.molarMass();
    node["Electronegativity"] = obj.electronegativity();
    node["Tags"]              = obj.tags();
}

auto operator>>(const yaml& node, Element& obj) -> void
{
    Element::Args args;
    node.at("Symbol").to(args.symbol);
    node.at("Name").to(args.name);
    node.at("AtomicNumber").to(args.atomic_number);
    node.at("AtomicWeight").to(args.atomic_weight);
    node.at("Electronegativity").to(args.electronegativity);
    node.at("Tags").to(args.tags);
    obj = Element(args);
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

}

auto operator>>(const yaml& node, Species& obj) -> void
{

}

} // namespace YAML
