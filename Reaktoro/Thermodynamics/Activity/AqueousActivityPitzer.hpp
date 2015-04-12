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

// C++ includes
#include <string>

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Activity/AqueousActivity.hpp>

namespace Reaktoro {

/// Create the aqueous activity function of the solvent species H<sub>2</sub>O(l) based on the Pitzer model
///
/// **References**
///   1. Harvie, C.E., Møller, N., Weare, J.H. (1984). The prediction of mineral solubilities in natural waters:
///      The Na-K-Mg-Ca-H-Cl-SO4-OH-HCO3-CO3-CO2-H2O system to high ionic strengths at 25°C.
///      Geochimica et Cosmochimica Acta, 48(4), 723–751.
///   2. Harvie, C.E., Weare, J.H. (1980). The prediction of mineral soluhilities in natural waters: the system
///      from zero to high concentration at 25°C. Geochimica et Cosmochimica Acta, 44(7), 981–997.
///   3. Pitzer, K. S. (1975). Thermodynamics of electrolytes. V. effects of higher-order electrostatic terms.
///      Journal of Solution Chemistry, 4(3), 249–265.
///
/// @param mixture The aqueous mixture instance containing the solvent species H<sub>2</sub>O(l)
/// @return The aqueous activity function of species H<sub>2</sub>O(l)
/// @see AqueousMixture, AqueousActivityFunction
auto aqueousActivityPitzerWater(const AqueousMixture& mixture) -> AqueousActivityFunction;

/**
 * Creates the aqueous activity function of a charged aqueous species based on the Pitzer model
 *
 * **References**
 *   1. Harvie, C.E., Møller, N., Weare, J.H. (1984). The prediction of mineral solubilities in natural waters:
 *      The Na-K-Mg-Ca-H-Cl-SO4-OH-HCO3-CO3-CO2-H2O system to high ionic strengths at 25°C.
 *      Geochimica et Cosmochimica Acta, 48(4), 723–751.
 *   2. Harvie, C.E., Weare, J.H. (1980). The prediction of mineral soluhilities in natural waters: the system
 *      from zero to high concentration at 25°C. Geochimica et Cosmochimica Acta, 44(7), 981–997.
 *   3. Pitzer, K. S. (1975). Thermodynamics of electrolytes. V. effects of higher-order electrostatic terms.
 *      Journal of Solution Chemistry, 4(3), 249–265.
 *
 * @param species The name of the aqueous species
 * @param mixture The aqueous mixture instance containing the aqueous species
 * @return The aqueous activity function of the aqueous species
 * @see AqueousMixture, AqueousActivity
 */
auto aqueousActivityPitzerCharged(const std::string& species, const AqueousMixture& mixture) -> AqueousActivityFunction;

/**
 * Creates the aqueous activity function of a neutral species based on the Pitzer model
 *
 * **References**
 *   1. Harvie, C.E., Møller, N., Weare, J.H. (1984). The prediction of mineral solubilities in natural waters:
 *      The Na-K-Mg-Ca-H-Cl-SO4-OH-HCO3-CO3-CO2-H2O system to high ionic strengths at 25°C.
 *      Geochimica et Cosmochimica Acta, 48(4), 723–751.
 *   2. Harvie, C.E., Weare, J.H. (1980). The prediction of mineral soluhilities in natural waters: the system
 *      from zero to high concentration at 25°C. Geochimica et Cosmochimica Acta, 44(7), 981–997.
 *   3. Pitzer, K. S. (1975). Thermodynamics of electrolytes. V. effects of higher-order electrostatic terms.
 *      Journal of Solution Chemistry, 4(3), 249–265.
 *
 * @param species The name of the aqueous species
 * @param mixture The aqueous mixture instance containing the aqueous species
 * @return The aqueous activity function of aqueous species
 * @see AqueousMixture, AqueousActivity
 */
auto aqueousActivityPitzerNeutral(const std::string& species, const AqueousMixture& mixture) -> AqueousActivityFunction;

} // namespace Reaktoro
