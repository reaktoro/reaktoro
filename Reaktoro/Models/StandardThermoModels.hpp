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

#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelConstant.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelHKF.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelHollandPowell.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelInterpolation.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelMaierKelley.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelMineralHKF.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelNasa.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelWaterHKF.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelYAML.hpp>

#include <Reaktoro/Models/StandardThermoModels/Support/SpeciesElectroProps.hpp>
#include <Reaktoro/Models/StandardThermoModels/Support/SpeciesElectroPropsHKF.hpp>

/// @defgroup StandardThermoModels Standard Thermodynamic Models
/// @ingroup Models
/// @brief The module in Reaktoro in which standard thermodynamic models for species and reactions are implemented.
