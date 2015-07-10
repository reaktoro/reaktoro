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
#include <complex>
#include <functional>
#include <tuple>

namespace Reaktoro {

/// Define a type that describes the roots of a cubic equation
using CubicRoots = std::tuple<std::complex<double>, std::complex<double>, std::complex<double>>;

/// Calculate the roots of a cubic equation using Cardano's method.
/// The calculation uses the approach presented in:
/// *Nickalls, R. W. D. (2012). A New Approach to Solving
/// the Cubic: Cardan’s Solution Revealed. The Mathematical
/// Gazette, 77(480), 354–359* to solve the cubic equation:
/// \f$ ax^{3}+bx^{2}+cx+d=0 \f$.
/// @param a The coefficient @c a of the cubic equation
/// @param b The coefficient @c b of the cubic equation
/// @param c The coefficient @c c of the cubic equation
/// @param d The coefficient @c d of the cubic equation
/// @return The three roots \f$ (r_1, r_2, r_3) \f$ that
/// solve the cubic equation, where \f$ |r_1| \geq |r_2| \geq |r_3| \f$.
auto cardano(double a, double b, double c, double d) -> CubicRoots;

/// Calculate the root of a non-linear function using Newton's method.
/// @param f The non-linear function \f$ f(x) \f$.
/// @param f The derivative \f$ f^{\prime}(x) \f$ of the non-linear function \f$ f(x) \f$.
/// @param x0 The initial guess for the iterative root calculation.
/// @param epsilon The tolerance used in \f$ |f(x)| < \epsilon \f$ to check convergence.
/// @param maxiter The maximum number of iterations.
auto newton(const std::function<double(double)>& f,
            const std::function<double(double)>& dfdx,
            double x0, double epsilon, unsigned maxiter) -> double;

} // namespace Reaktoro
