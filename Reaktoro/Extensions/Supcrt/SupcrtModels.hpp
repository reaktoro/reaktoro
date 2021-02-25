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
#include <Reaktoro/Extensions/Supcrt/SupcrtParams.hpp>

namespace Reaktoro {

// Forward declarations
struct SpeciesElectroState;
struct SpeciesThermoState;
struct WaterElectroState;
struct WaterThermoState;

/// Calculate the thermodynamic state of solvent water using the HKF model.
auto supcrtStandardThermoPropsSolventHKF(real T, real P, const SupcrtParamsAqueousSolventHKF& params, const WaterThermoState& wts) -> SpeciesThermoState;

/// Calculate the thermodynamic state of an aqueous solute using the HKF model.
auto supcrtStandardThermoPropsSoluteHKF(real T, real P, const SupcrtParamsAqueousSoluteHKF& params, const SpeciesElectroState& aes, const WaterElectroState& wes) -> SpeciesThermoState;

/// Calculate the thermodynamic state of a fluid species using the Maier-Kelly model.
auto supcrtStandardThermoPropsMaierKelly(real T, real P, const SupcrtParamsMaierKelly& params) -> SpeciesThermoState;

/// Calculate the thermodynamic state of a mineral species using the Maier-Kelly-HKF model.
auto supcrtStandardThermoPropsMaierKellyHKF(real T, real P, const SupcrtParamsMaierKellyHKF& params) -> SpeciesThermoState;

} // namespace Reaktoro
