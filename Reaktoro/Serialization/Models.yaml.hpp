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

// Forward declarations
class yaml;

struct ReactionThermoModelParamsConstLgK;
struct ReactionThermoModelParamsGemsLgK;
struct ReactionThermoModelParamsPhreeqcLgK;
struct ReactionThermoModelParamsVantHoff;
struct StandardThermoModelParamsHKF;
struct StandardThermoModelParamsHollandPowell;
struct StandardThermoModelParamsMaierKelley;
struct StandardThermoModelParamsMineralHKF;
struct StandardThermoModelParamsWaterHKF;

auto operator<<(yaml& node, const ReactionThermoModelParamsConstLgK& obj) -> void;
auto operator>>(const yaml& node, ReactionThermoModelParamsConstLgK& obj) -> void;

auto operator<<(yaml& node, const ReactionThermoModelParamsGemsLgK& obj) -> void;
auto operator>>(const yaml& node, ReactionThermoModelParamsGemsLgK& obj) -> void;

auto operator<<(yaml& node, const ReactionThermoModelParamsPhreeqcLgK& obj) -> void;
auto operator>>(const yaml& node, ReactionThermoModelParamsPhreeqcLgK& obj) -> void;

auto operator<<(yaml& node, const ReactionThermoModelParamsVantHoff& obj) -> void;
auto operator>>(const yaml& node, ReactionThermoModelParamsVantHoff& obj) -> void;

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
