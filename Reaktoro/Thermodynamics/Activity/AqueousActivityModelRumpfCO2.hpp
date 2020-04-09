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
#include <Reaktoro/Thermodynamics/Activity/AqueousActivityModel.hpp>

namespace Reaktoro {

/// Create the function for the calculation of ln activity coefficient of CO<sub>2</sub>(aq), based on the model of Rumpf et al. (1994).
/// The model is documented in the paper *Rumpf, B., Nicolaisen, H., Ocal, C., &
/// Maurer, G. (1994). Solubility of carbon dioxide in aqueous mixtures of sodium
/// chloride: Experimental results and correlation. Journal of Solution Chemistry,
// 23(3), 431–448*.
/// @param mixture The aqueous mixture instance
/// @see AqueousMixture, AqueousActivityModel
auto aqueousActivityModelRumpfCO2(const AqueousMixture& mixture) -> AqueousActivityModel;

} // namespace Reaktoro
