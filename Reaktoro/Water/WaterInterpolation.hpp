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

namespace Reaktoro {

/// Compute the density of water at given a temperature and pressure using quadratic interpolation.
/// The interpolation is performed using collected water properties at temperatures and pressures
/// shown in Table 13.2 of *Wagner, W., Pruss, A. (2002). The IAPWS Formulation 1995 for the
/// Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use. Journal of
/// Physical and Chemical Reference Data, 31(2), 387. https://doi.org/10.1063/1.1461829*.
/// @param T The temperature value (in K)
/// @param P The pressure value (in Pa)
auto waterDensityWagnerPrussInterp(real const& T, real const& P) -> real;

} // namespace Reaktoro
