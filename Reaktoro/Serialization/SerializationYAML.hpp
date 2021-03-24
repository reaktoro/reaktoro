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

namespace Reaktoro {

// Forward declarations (Core)
class ChemicalFormula;
class ChemicalSystem;
class Element;
class Param;
class Params;
class Phase;
class Species;

// Forward declarations (StandardThermoModel)
struct StandardThermoModelParamsHKF;
struct StandardThermoModelParamsHollandPowell;
struct StandardThermoModelParamsMaierKelley;
struct StandardThermoModelParamsMineralHKF;
struct StandardThermoModelParamsWaterHKF;

} // namespace Reaktoro

namespace YAML {

using namespace Reaktoro;

// Forward declaration
class Node;

//=====================================================================================================================
// Common
//=====================================================================================================================

auto operator<<(Node& node, const real& obj) -> void;
auto operator>>(const Node& node, real& obj) -> void;

//=====================================================================================================================
// Core
//=====================================================================================================================

auto operator<<(Node& node, const ChemicalFormula& obj) -> void;
auto operator>>(const Node& node, ChemicalFormula& obj) -> void;

auto operator<<(Node& node, const ChemicalSystem& obj) -> void;
auto operator>>(const Node& node, ChemicalSystem& obj) -> void;

auto operator<<(Node& node, const Element& obj) -> void;
auto operator>>(const Node& node, Element& obj) -> void;

auto operator<<(Node& node, const Param& obj) -> void;
auto operator>>(const Node& node, Param& obj) -> void;

auto operator<<(Node& node, const Params& obj) -> void;
auto operator>>(const Node& node, Params& obj) -> void;

auto operator<<(Node& node, const Phase& obj) -> void;
auto operator>>(const Node& node, Phase& obj) -> void;

auto operator<<(Node& node, const Species& obj) -> void;
auto operator>>(const Node& node, Species& obj) -> void;

//=====================================================================================================================
// StandardThermoModel
//=====================================================================================================================

auto operator<<(Node& node, const StandardThermoModelParamsHKF& obj) -> void;
auto operator>>(const Node& node, StandardThermoModelParamsHKF& obj) -> void;

auto operator<<(Node& node, const StandardThermoModelParamsHollandPowell& obj) -> void;
auto operator>>(const Node& node, StandardThermoModelParamsHollandPowell& obj) -> void;

auto operator<<(Node& node, const StandardThermoModelParamsMaierKelley& obj) -> void;
auto operator>>(const Node& node, StandardThermoModelParamsMaierKelley& obj) -> void;

auto operator<<(Node& node, const StandardThermoModelParamsMineralHKF& obj) -> void;
auto operator>>(const Node& node, StandardThermoModelParamsMineralHKF& obj) -> void;

auto operator<<(Node& node, const StandardThermoModelParamsWaterHKF& obj) -> void;
auto operator>>(const Node& node, StandardThermoModelParamsWaterHKF& obj) -> void;

} // namespace YAML
