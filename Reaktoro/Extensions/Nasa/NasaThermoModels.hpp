// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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
#include <Reaktoro/Core/StandardThermoProps.hpp>
#include <Reaktoro/Extensions/Nasa/NasaSpecies.hpp>
#include <Reaktoro/Extensions/Nasa/NasaThermoData.hpp>

namespace Reaktoro {
namespace NasaUtils {

/// Return true if given thermodynamic parameters have continuous temperature intervals.
auto areTemperatureIntervalsContinuous(const Vec<NasaThermoData>& data) -> bool;

/// Return the minimum supported temperature (in K) for provided thermodynamic data.
auto minSupportedTemperature(const Vec<NasaThermoData>& data) -> real;

/// Return the maximum supported temperature (in K) for provided thermodynamic data.
auto maxSupportedTemperature(const Vec<NasaThermoData>& data) -> real;

/// Return the NasaThermoData object that is applicable for given temperature value (in K).
auto getNasaThermoDataForGivenTemperature(const Vec<NasaThermoData>& data, const real& T) -> NasaThermoData;

/// Compute the standard thermodynamic properties of a species with given NASA thermodynamic parameters.
auto computeStandardThermoProps(const NasaThermoData& params, const real& T) -> StandardThermoProps;

/// Compute the standard thermodynamic properties of a species with given NASA thermodynamic parameters.
auto computeStandardThermoProps(const NasaSpecies& species, const real& T) -> StandardThermoProps;

} // namespace NasaUtils

/// Return a standard thermodynamic model for a given chemical species with NASA thermodynamic parameters.
auto StandardThermoModelNasa(const NasaSpecies& species) -> StandardThermoModel;

} // namespace Reaktoro
