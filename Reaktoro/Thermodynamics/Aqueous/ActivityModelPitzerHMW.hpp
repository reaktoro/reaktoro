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
#include <Reaktoro/Core/ActivityModel.hpp>

namespace Reaktoro {

/// Return the activity model for aqueous phases based on Harvie-Møller-Weare Pitzer's formulation.
/// **References:**
///   - Harvie, C.E., Møller, N., Weare, J.H. (1984). The prediction of mineral
///     solubilities in natural waters: The
///     Na-K-Mg-Ca-H-Cl-SO4-OH-HCO3-CO3-CO2-H2O system to high ionic strengths
///     at 25°C. Geochimica et Cosmochimica Acta, 48(4), 723–751.
///   - Harvie, C.E., Weare, J.H. (1980). The prediction of mineral
///     soluhilities in natural waters: the system from zero to high
///     concentration at 25°C. Geochimica et Cosmochimica Acta, 44(7), 981–997.
///   - Pitzer, K. S. (1975). Thermodynamics of electrolytes. V. effects of
///     higher-order electrostatic terms. Journal of Solution Chemistry, 4(3),
///     249–265.
///   - Helgeson, H. C., Kirkham, D. H., Flowers, G. C. (1981). Theoretical
///     prediction of the thermodynamic behavior of aqueous electrolytes at
///     high pressures and temperatures: IV. Calculation of activity
///     coefficients, osmotic coefficients, and apparent molal and standard and
///     relative partial molal properties to 600°C. American Journal of
///     Science, 281(10), 1249–1516.
/// @ingroup ActivityModels
auto ActivityModelPitzerHMW() -> ActivityModelGenerator;

} // namespace Reaktoro
