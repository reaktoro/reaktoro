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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

// Reaktor includes
#include <Reaktor/Common/VectorResult.hpp>

namespace Reaktor {

/**
 * Defines the functional signature of a vector-valued activity function
 * @param T The temperature for the activity calculation (in units of K)
 * @param P The pressure for the activity calculation (in units of Pa)
 * @param n The molar composition of the species (in units of mol)
 * @return The activities of the species and their molar derivatives
 * @see Phase, VectorResult
 * @ingroup Core
 */
using ActivityFn = std::function<VectorResult(double T, double P, const Vector& n)>;

/**
 * Defines the functional signature of a concentration function
 * @param n The molar composition of the species (in units of mol)
 * @return The concentration of the species
 * @see Phase, Vector
 * @ingroup Core
 */
using ConcentrationFn = std::function<Vector(const Vector& n)>;

/**
 * Defines the function signature of the standart chemical potential of a species
 * @param T The temperature for the function evaluation (in units of K)
 * @param P The pressure for the function evaluation (in units of Pa)
 * @return The standart chemical potential of a species
 * @ingroup Common
 */
using ChemicalPotentialFn = std::function<double(double, double)>;

/**
 * Defines the function signature of the density of a species (in units of kg/m3)
 * @param T The temperature for the function evaluation (in units of K)
 * @param P The pressure for the function evaluation (in units of Pa)
 * @return The density of a species (in units of kg/m3)
 * @ingroup Common
 */
using DensityFn = std::function<double(double, double)>;

/**
 * Defines the function signature of an equilibrium constant of a reaction
 * @param T The temperature for the calculation (in units of K)
 * @param P The pressure for the calculation (in units of Pa)
 * @return The equilibrium constant of the reaction
 */
using EquilibriumConstantFn = std::function<double(double, double)>;

/**
 * Defines the functional signature of a kinetic rate function
 * @param T The temperature for the rate calculation (in units of K)
 * @param P The pressure for the rate calculation (in units of Pa)
 * @param n The molar abundance of the species (in units of mol)
 * @param a The activities of the species and their molar derivatives
 * @return The kinetic rate of the reaction and its molar derivatives
 * @see Reaction, ScalarResult, VectorResult
 * @ingroup Core
 */
using ReactionRateFn = std::function<ScalarResult(double T, double P, const Vector& n, const VectorResult& a)>;

} // namespace Reaktor


