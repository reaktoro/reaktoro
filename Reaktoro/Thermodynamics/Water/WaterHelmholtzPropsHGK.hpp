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
#include <Reaktoro/Common/Real.hpp>

namespace Reaktoro {

// Forward declarations
struct WaterHelmholtzProps;

/// Calculate the Helmholtz free energy state of water using the Haar--Gallagher--Kell (1984) equation of state
/// @param T The temperature of water (in units of K)
/// @param D The density of water (in units of kg/m3)
/// @return The Helmholtz free energy state of water
/// @see WaterHelmholtzProps
auto waterHelmholtzPropsHGK(real T, real D) -> WaterHelmholtzProps;

} // namespace Reaktoro
