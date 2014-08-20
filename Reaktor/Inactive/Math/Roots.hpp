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
#include <complex>
#include <tuple>

namespace Reaktor {

/**
 * Type that describes the roots of a cubic equation
 */
using CubicRoots = std::tuple<std::complex<double>, std::complex<double>, std::complex<double>>;

/**
 * Calculates the roots of a cubic polynomial
 *
 * The calculation uses the approach presented in:
 *  Nickalls, R. W. D. (2012). A New Approach to Solving
 *  the Cubic: Cardan’s Solution Revealed. The Mathematical
 *  Gazette, 77(480), 354–359
 * to solve the cubic equation: \f$ ax^{3}+bx^{2}+cx+d=0 \f$
 *
 * @param a The coefficient @c a of the cubic equation
 * @param b The coefficient @c b of the cubic equation
 * @param c The coefficient @c c of the cubic equation
 * @param d The coefficient @c d of the cubic equation
 * @return The three roots that solve the cubic equation
 */
auto cubicRoots(double a, double b, double c, double d) -> CubicRoots;

} /* namespace Reaktor */
