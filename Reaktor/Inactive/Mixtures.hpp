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

// Reaktor includes
#include <Reaktor/Mixtures/AqueousMixture.hpp>
#include <Reaktor/Mixtures/GaseousMixture.hpp>
#include <Reaktor/Mixtures/GeneralMixture.hpp>
#include <Reaktor/Mixtures/MineralMixture.hpp>

/**
 * @defgroup Mixtures Mixtures
 *
 * Provides the definitions of mixtures of species
 *
 * A mixture is defined as a collection of species. For example, an AqueousMixture is
 * a collection of AqueousSpecies objects. Each mixture should provide the necessary
 * operations for activity calculations. Thus, a AqueousMixture implements, for instance,
 * methods for the calculation of molar fractions, molalities, effective and stoichiometric
 * ionic strengths, and so forth.
 */
