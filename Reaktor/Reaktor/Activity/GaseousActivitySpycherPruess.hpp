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
#include <Reaktor/Activity/GaseousActivity.hpp>

namespace Reaktor {

/// Create the gaseous activity functions of species H<sub>2</sub>O(g) and CO<sub>2</sub>(g) based on the model of Spycher and Pruess (2003)
///
/// @b References
/// 1. Spycher, N., Pruess, K., Ennis-King, J. (2003). CO2--H2O mixtures in the geological sequestration of CO2. I. Assessment and calculation of mutual solubilities from 12 to 100°C and up to 600 bar. Geochimica et Cosmochimica Acta, 67(16), 3015–3031. doi:10.1016/S0016-7037(03)00273-4
///
/// @param mixture The gaseous mixture instance
/// @return The gaseous activity functions of species H<sub>2</sub>O(g) and CO<sub>2</sub>(g) (in this order)
/// @see GaseousMixture, GaseousActivity
auto gaseousActivitySpycherPruessH2OCO2(const GaseousMixture& mixture) -> std::vector<GaseousActivity>;

} // namespace Reaktor
