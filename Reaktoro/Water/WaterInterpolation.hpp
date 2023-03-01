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
#include <Reaktoro/Core/StateOfMatter.hpp>

namespace Reaktoro {

// Forward declarations
struct WaterThermoProps;

/// Compute the density of water (in kg/m3) at given a temperature and pressure using quadratic interpolation.
/// The interpolation is performed using pre-computed water properties at temperatures and pressures
/// shown in Table 13.2 of *Wagner, W., Pruss, A. (2002). The IAPWS Formulation 1995 for the
/// Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use. Journal of
/// Physical and Chemical Reference Data, 31(2), 387. https://doi.org/10.1063/1.1461829*.
/// @param T The temperature value (in K)
/// @param P The pressure value (in Pa)
/// @param som The desired state of matter for water (the actual state of matter may end up being different!)
auto waterDensityWagnerPrussInterp(real const& T, real const& P, StateOfMatter som) -> real;

/// Compute the thermodynamic properties of water at given a temperature and pressure using quadratic interpolation.
/// The interpolation is performed using pre-computed water properties at temperatures and pressures
/// shown in Table 13.2 of *Wagner, W., Pruss, A. (2002). The IAPWS Formulation 1995 for the
/// Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use. Journal of
/// Physical and Chemical Reference Data, 31(2), 387. https://doi.org/10.1063/1.1461829*.
/// @param T The temperature value (in K)
/// @param P The pressure value (in Pa)
/// @param som The desired state of matter for water (the actual state of matter may end up being different!)
auto waterThermoPropsWagnerPrussInterp(real const& T, real const& P, StateOfMatter som) -> WaterThermoProps;

/// Return the pre-computed thermodynamic properties of water using Wagner and Pruss (2002) equation of state used for interpolation.
/// The interpolation data consists of pre-computed water properties at temperatures and pressures
/// shown in Table 13.2 of *Wagner, W., Pruss, A. (2002). The IAPWS Formulation 1995 for the
/// Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use. Journal of
/// Physical and Chemical Reference Data, 31(2), 387. https://doi.org/10.1063/1.1461829*.
/// @note The first call to this function will trigger a parsing operation of embedded text data.
/// @param som The desired state of matter for water (the actual state of matter may end up being different!)
auto waterThermoPropsWagnerPrussInterpData(StateOfMatter som) -> Vec<Vec<WaterThermoProps>> const&;

} // namespace Reaktoro
