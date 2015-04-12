// Reaktoro is a C++ library for computational reaction modelling.
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

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Activity/AqueousActivity.hpp>

namespace Reaktoro {

/// Create the aqueous activity function of species CO<sub>2</sub>(aq) based on the model of Duan and Sun (2003)
///
/// @b References
/// 1. Duan, Z., Sun, R. (2003). An improved model calculating CO2 solubility in pure water and aqueous NaCl mixtures from 273 to 533 K and from 0 to 2000 bar. Chemical Geology, 193(3-4), 257–271. doi:10.1016/S0009-2541(02)00263-2
///
/// @param mixture The aqueous mixture instance
/// @return The aqueous activity function of species CO<sub>2</sub>(aq)
/// @see AqueousMixture, AqueousActivityFunction
auto aqueousActivityDuanSunCO2(const AqueousMixture& mixture) -> AqueousActivityFunction;

} // namespace Reaktoro
