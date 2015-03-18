// Reaktor is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#pragma once

// Reaktor includes
#include <Reaktor/Thermodynamics/Activity/GaseousActivity.hpp>

namespace Reaktor {

/// Create the gaseous activity functions of species H<sub>2</sub>O(g), CO<sub>2</sub>(g) and CH<sub>4</sub>(g) based on the model of Spycher and Pruess (2003)
///
/// @b References
/// 1. Spycher, N., Reed, M. (1988). Fugacity coefficients of H2, CO2, CH4, H2O and of H2O--CO2--CH4 solutions: A virial equation treatment for moderate pressures and temperatures applicable to calculations of hydrothermal boiling. Geochimica et Cosmochimica Acta, 52(3), 739–749. doi:10.1016/0016-7037(88)90334-1
///
/// @param solution The gaseous solution instance
/// @return The gaseous activity functions of species H<sub>2</sub>O(g), CO<sub>2</sub>(g) and CH<sub>4</sub>(g) (in this order)
/// @see GaseousSolution, GaseousActivity
auto gaseousActivitySpycherReedH2OCO2CH4(const GaseousSolution& solution) -> std::vector<GaseousActivity>;

} // namespace Reaktor
