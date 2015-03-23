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
#include <Reaktor/Thermodynamics/Activity/AqueousActivity.hpp>

namespace Reaktor {

/// Create the aqueous activity function of species CO<sub>2</sub>(aq) based on the model of Rumpf et al. (1994)
///
/// @b References
/// 1. Rumpf, B., Nicolaisen, H., Ocal, C., & Maurer, G. (1994). Solubility of carbon dioxide in aqueous mixtures of sodium chloride: Experimental results and correlation. Journal of Solution Chemistry, 23(3), 431–448. doi:10.1007/BF00973113
///
/// @param mixture The aqueous mixture instance
/// @return The aqueous activity function of species CO<sub>2</sub>(aq)
/// @see AqueousMixture, AqueousActivityFunction
auto aqueousActivityRumpfCO2(const AqueousMixture& mixture) -> AqueousActivityFunction;

} // namespace Reaktor
