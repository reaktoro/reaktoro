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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/AggregateState.hpp>

namespace Reaktoro {

// Forward declarations
class yaml;
class ChemicalFormula;
class ChemicalSystem;
class Element;
class ElementalComposition;
class FormationReaction;
class Param;
class Params;
class Phase;
class Species;
template<typename Signature> class Model;

auto operator<<(yaml& node, const AggregateState& obj) -> void;
auto operator>>(const yaml& node, AggregateState& obj) -> void;

auto operator<<(yaml& node, const ChemicalFormula& obj) -> void;
auto operator>>(const yaml& node, ChemicalFormula& obj) -> void;

auto operator<<(yaml& node, const ChemicalSystem& obj) -> void;
auto operator>>(const yaml& node, ChemicalSystem& obj) -> void;

auto operator<<(yaml& node, const Element& obj) -> void;
auto operator>>(const yaml& node, Element& obj) -> void;

auto operator<<(yaml& node, const ElementalComposition& obj) -> void;
auto operator>>(const yaml& node, ElementalComposition& obj) -> void;

auto operator<<(yaml& node, const FormationReaction& obj) -> void;
auto operator>>(const yaml& node, FormationReaction& obj) -> void;

auto operator<<(yaml& node, const Param& obj) -> void;
auto operator>>(const yaml& node, Param& obj) -> void;

auto operator<<(yaml& node, const Params& obj) -> void;
auto operator>>(const yaml& node, Params& obj) -> void;

auto operator<<(yaml& node, const Phase& obj) -> void;
auto operator>>(const yaml& node, Phase& obj) -> void;

auto operator<<(yaml& node, const Species& obj) -> void;
auto operator>>(const yaml& node, Species& obj) -> void;

} // namespace Reaktoro

// REAKTORO_DECLARE_CONVERT_YAML(ChemicalFormula);

// namespace YAML {

// template<>
// struct convert<ChemicalFormula>
// {
//     static auto encode(const ChemicalFormula& obj);
//     static auto decode(const Node& node, ChemicalFormula& obj);
// };

// } // namespace YAML
