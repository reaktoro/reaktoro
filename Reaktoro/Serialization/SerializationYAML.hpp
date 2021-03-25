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

// Forward declarations (Common)
class yaml;

// Forward declarations (Core)
class ChemicalFormula;
class ChemicalSystem;
class Element;
class Param;
class Params;
class Phase;
class Species;

template<typename Signature>
class Model;

// Forward declarations (StandardThermoModel)
struct StandardThermoModelParamsHKF;
struct StandardThermoModelParamsHollandPowell;
struct StandardThermoModelParamsMaierKelley;
struct StandardThermoModelParamsMineralHKF;
struct StandardThermoModelParamsWaterHKF;

} // namespace Reaktoro

namespace Reaktoro {

//=====================================================================================================================
// Common
//=====================================================================================================================

auto operator<<(yaml& node, const real& obj) -> void;
auto operator>>(const yaml& node, real& obj) -> void;

//=====================================================================================================================
// Core
//=====================================================================================================================

auto operator<<(yaml& node, const ChemicalFormula& obj) -> void;
auto operator>>(const yaml& node, ChemicalFormula& obj) -> void;

auto operator<<(yaml& node, const ChemicalSystem& obj) -> void;
auto operator>>(const yaml& node, ChemicalSystem& obj) -> void;

auto operator<<(yaml& node, const Element& obj) -> void;
auto operator>>(const yaml& node, Element& obj) -> void;

auto operator<<(yaml& node, const Param& obj) -> void;
auto operator>>(const yaml& node, Param& obj) -> void;

auto operator<<(yaml& node, const Params& obj) -> void;
auto operator>>(const yaml& node, Params& obj) -> void;

auto operator<<(yaml& node, const Phase& obj) -> void;
auto operator>>(const yaml& node, Phase& obj) -> void;

auto operator<<(yaml& node, const Species& obj) -> void;
auto operator>>(const yaml& node, Species& obj) -> void;

// template<typename Signature>
// auto operator<<(yaml& node, const Model<Signature>& obj) -> void
// {
//     node = obj.yaml();
// }


// auto operator<<(yaml& node, const std::any& params) -> void

// // template<typename Signature>
// // auto operator>>(const yaml& node, Model<Signature>& obj) -> void
// // {

// // }

//=====================================================================================================================
// StandardThermoModel
//=====================================================================================================================

auto operator<<(yaml& node, const StandardThermoModelParamsHKF& obj) -> void;
auto operator>>(const yaml& node, StandardThermoModelParamsHKF& obj) -> void;

auto operator<<(yaml& node, const StandardThermoModelParamsHollandPowell& obj) -> void;
auto operator>>(const yaml& node, StandardThermoModelParamsHollandPowell& obj) -> void;

auto operator<<(yaml& node, const StandardThermoModelParamsMaierKelley& obj) -> void;
auto operator>>(const yaml& node, StandardThermoModelParamsMaierKelley& obj) -> void;

auto operator<<(yaml& node, const StandardThermoModelParamsMineralHKF& obj) -> void;
auto operator>>(const yaml& node, StandardThermoModelParamsMineralHKF& obj) -> void;

auto operator<<(yaml& node, const StandardThermoModelParamsWaterHKF& obj) -> void;
auto operator>>(const yaml& node, StandardThermoModelParamsWaterHKF& obj) -> void;

} // namespace YAML
