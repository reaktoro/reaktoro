// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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
#include <Reaktoro/Common/YAML.hpp>

namespace Reaktoro {

//======================================================================
// ReactionThermoModel Types
//======================================================================
struct ReactionThermoModelParamsConstLgK;
struct ReactionThermoModelParamsGemsLgK;
struct ReactionThermoModelParamsPhreeqcLgK;
struct ReactionThermoModelParamsVantHoff;

REAKTORO_YAML_ENCODE_DECLARE(ReactionThermoModelParamsConstLgK);
REAKTORO_YAML_DECODE_DECLARE(ReactionThermoModelParamsConstLgK);

REAKTORO_YAML_ENCODE_DECLARE(ReactionThermoModelParamsGemsLgK);
REAKTORO_YAML_DECODE_DECLARE(ReactionThermoModelParamsGemsLgK);

REAKTORO_YAML_ENCODE_DECLARE(ReactionThermoModelParamsPhreeqcLgK);
REAKTORO_YAML_DECODE_DECLARE(ReactionThermoModelParamsPhreeqcLgK);

REAKTORO_YAML_ENCODE_DECLARE(ReactionThermoModelParamsVantHoff);
REAKTORO_YAML_DECODE_DECLARE(ReactionThermoModelParamsVantHoff);

//======================================================================
// StandardThermoModel Types
//======================================================================
struct StandardThermoModelParamsConstant;
struct StandardThermoModelParamsHKF;
struct StandardThermoModelParamsHollandPowell;
struct StandardThermoModelParamsMaierKelley;
struct StandardThermoModelParamsMineralHKF;
struct StandardThermoModelParamsWaterHKF;

REAKTORO_YAML_ENCODE_DECLARE(StandardThermoModelParamsConstant);
REAKTORO_YAML_DECODE_DECLARE(StandardThermoModelParamsConstant);

REAKTORO_YAML_ENCODE_DECLARE(StandardThermoModelParamsHKF);
REAKTORO_YAML_DECODE_DECLARE(StandardThermoModelParamsHKF);

REAKTORO_YAML_ENCODE_DECLARE(StandardThermoModelParamsHollandPowell);
REAKTORO_YAML_DECODE_DECLARE(StandardThermoModelParamsHollandPowell);

REAKTORO_YAML_ENCODE_DECLARE(StandardThermoModelParamsMaierKelley);
REAKTORO_YAML_DECODE_DECLARE(StandardThermoModelParamsMaierKelley);

REAKTORO_YAML_ENCODE_DECLARE(StandardThermoModelParamsMineralHKF);
REAKTORO_YAML_DECODE_DECLARE(StandardThermoModelParamsMineralHKF);

REAKTORO_YAML_ENCODE_DECLARE(StandardThermoModelParamsWaterHKF);
REAKTORO_YAML_DECODE_DECLARE(StandardThermoModelParamsWaterHKF);

//======================================================================
// StandardVolumeModel Types
//======================================================================
struct StandardVolumeModelParamsConstant;

REAKTORO_YAML_ENCODE_DECLARE(StandardVolumeModelParamsConstant);
REAKTORO_YAML_DECODE_DECLARE(StandardVolumeModelParamsConstant);

} // namespace YAML
