/*
 * Reaktor is a C++ library for computational reaction modelling.
 *
 * Copyright (C) 2014 Allan Leal
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

// C++ includes
#include <string>

// Reaktor includes
#include <Reaktor/Activity/AqueousActivity.hpp>

namespace Reaktor {

// Forward declarations
class AqueousMixture;

/**
 * Creates the aqueous activity function of a charged species based on the HKF model
 *
 * **References**
 *   1. Helgeson, H. C., Kirkham, D. H., Flowers, G. C. (1981). Theoretical prediction of the thermodynamic behavior
 *      of aqueous electrolytes at high pressures and temperatures: IV. Calculation of activity coefficients, osmotic
 *      coefficients, and apparent molal and standard and relative partial molal properties to 600°C.
 *      American Journal of Science, 281(10), 1249–1516.
 *
 * @param species The name of the charged species
 * @param mixture The aqueous mixture instance containing the charged species
 * @return The aqueous activity function of the charged species
 * @see AqueousMixture, AqueousActivity
 */
auto aqueousActivityHKFCharged(const std::string& species, const AqueousMixture& mixture) -> AqueousActivity;

/**
 * Creates the aqueous activity function of the solvent species H<sub>2</sub>O(l) based on the HKF model
 *
 * **References**
 *   1. Helgeson, H. C., Kirkham, D. H., Flowers, G. C. (1981). Theoretical prediction of the thermodynamic behavior
 *      of aqueous electrolytes at high pressures and temperatures: IV. Calculation of activity coefficients, osmotic
 *      coefficients, and apparent molal and standard and relative partial molal properties to 600°C.
 *      American Journal of Science, 281(10), 1249–1516.
 *
 * @param mixture The aqueous mixture instance containing the solvent species H<sub>2</sub>O(l)
 * @return The aqueous activity function of species H<sub>2</sub>O(l)
 * @see AqueousMixture, AqueousActivity
 */
auto aqueousActivityHKFWater(const AqueousMixture& mixture) -> AqueousActivity;

} // namespace Reaktor
