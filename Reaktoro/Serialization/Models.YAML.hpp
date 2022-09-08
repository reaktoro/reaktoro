// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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
// StandardThermoModelParams Types
//======================================================================
struct StandardThermoModelParamsConstant;
struct StandardThermoModelParamsHKF;
struct StandardThermoModelParamsHollandPowell;
struct StandardThermoModelParamsInterpolation;
struct StandardThermoModelParamsMaierKelley;
struct StandardThermoModelParamsMineralHKF;
struct StandardThermoModelParamsNasa;
struct StandardThermoModelParamsWaterHKF;

REAKTORO_YAML_ENCODE_DECLARE(StandardThermoModelParamsConstant);
REAKTORO_YAML_DECODE_DECLARE(StandardThermoModelParamsConstant);

REAKTORO_YAML_ENCODE_DECLARE(StandardThermoModelParamsHKF);
REAKTORO_YAML_DECODE_DECLARE(StandardThermoModelParamsHKF);

REAKTORO_YAML_ENCODE_DECLARE(StandardThermoModelParamsHollandPowell);
REAKTORO_YAML_DECODE_DECLARE(StandardThermoModelParamsHollandPowell);

REAKTORO_YAML_ENCODE_DECLARE(StandardThermoModelParamsInterpolation);
REAKTORO_YAML_DECODE_DECLARE(StandardThermoModelParamsInterpolation);

REAKTORO_YAML_ENCODE_DECLARE(StandardThermoModelParamsMaierKelley);
REAKTORO_YAML_DECODE_DECLARE(StandardThermoModelParamsMaierKelley);

REAKTORO_YAML_ENCODE_DECLARE(StandardThermoModelParamsMineralHKF);
REAKTORO_YAML_DECODE_DECLARE(StandardThermoModelParamsMineralHKF);

REAKTORO_YAML_ENCODE_DECLARE(StandardThermoModelParamsNasa);
REAKTORO_YAML_DECODE_DECLARE(StandardThermoModelParamsNasa);

REAKTORO_YAML_ENCODE_DECLARE(StandardThermoModelParamsWaterHKF);
REAKTORO_YAML_DECODE_DECLARE(StandardThermoModelParamsWaterHKF);

//======================================================================
// ReactionStandardThermoModelParams Types
//======================================================================
struct ReactionStandardThermoModelParamsConstLgK;
struct ReactionStandardThermoModelParamsGemsLgK;
struct ReactionStandardThermoModelParamsPhreeqcLgK;
struct ReactionStandardThermoModelParamsVantHoff;

REAKTORO_YAML_ENCODE_DECLARE(ReactionStandardThermoModelParamsConstLgK);
REAKTORO_YAML_DECODE_DECLARE(ReactionStandardThermoModelParamsConstLgK);

REAKTORO_YAML_ENCODE_DECLARE(ReactionStandardThermoModelParamsGemsLgK);
REAKTORO_YAML_DECODE_DECLARE(ReactionStandardThermoModelParamsGemsLgK);

REAKTORO_YAML_ENCODE_DECLARE(ReactionStandardThermoModelParamsPhreeqcLgK);
REAKTORO_YAML_DECODE_DECLARE(ReactionStandardThermoModelParamsPhreeqcLgK);

REAKTORO_YAML_ENCODE_DECLARE(ReactionStandardThermoModelParamsVantHoff);
REAKTORO_YAML_DECODE_DECLARE(ReactionStandardThermoModelParamsVantHoff);

//======================================================================
// StandardVolumeModelParams Types
//======================================================================
struct StandardVolumeModelParamsConstant;

REAKTORO_YAML_ENCODE_DECLARE(StandardVolumeModelParamsConstant);
REAKTORO_YAML_DECODE_DECLARE(StandardVolumeModelParamsConstant);

//======================================================================
// ReactionRateModelParams Types
//======================================================================
struct ReactionRateModelParamsPalandriKharaka;

REAKTORO_YAML_ENCODE_DECLARE(ReactionRateModelParamsPalandriKharaka);
REAKTORO_YAML_DECODE_DECLARE(ReactionRateModelParamsPalandriKharaka);

} // namespace YAML
