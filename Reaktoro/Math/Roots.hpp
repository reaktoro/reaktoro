// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

// C++ includes
#include <complex>
#include <functional>
#include <tuple>

namespace Reaktoro {

/// Define a type that describes the roots of a cubic equation
template<typename T>
using CubicRoots = std::tuple<std::complex<T>, std::complex<T>, std::complex<T>>;

/// Calculate the roots of a cubic equation using Cardano's method.
/// The calculation uses the approach presented in:
///
///    Nickalls, R. W. D. (2012). A New Approach to Solving the Cubic: Cardan’s
///    Solution Revealed. The Mathematical Gazette, 77(480), 354–359* to solve
///    the cubic equation: \f$ x^{3}+bx^{2}+cx+d=0 \f$.
///
/// @param b The coefficient *b* of the cubic equation
/// @param c The coefficient *c* of the cubic equation
/// @param d The coefficient *d* of the cubic equation
/// @return The three roots \f$ (r_1, r_2, r_3) \f$ that solve the cubic
/// equation, where \f$ |r_1| \geq |r_2| \geq |r_3| \f$. If they are all real,
/// they are returned in descending order (the first element is the highest).
template<typename T>
auto cardano(const T& b, const T& c, const T& d) -> CubicRoots<T>
{
    using std::abs;
    using std::acos;
    using std::cos;
    using std::pow;
    using std::sqrt;

    const T xn = -b/3.0;
    const T yn = xn*xn*xn + b*xn*xn + c*xn + d;

    const T d2    = (b*b - 3*c)/9.0;
    const T h2    = 4*d2*d2*d2;
    const T Delta = yn*yn - h2;

    std::complex<T> x1, x2, x3;

    if(Delta > 0.0)
    {
        const T sqrtDelta = sqrt(Delta);

        const T operand1 = (-yn + sqrtDelta) * 0.5;
        const T operand2 = (-yn - sqrtDelta) * 0.5;

        const T alpha = xn +
            pow(abs(operand1), 1.0/3) * abs(operand1)/operand1 +
            pow(abs(operand2), 1.0/3) * abs(operand2)/operand2;

        const T discr = b*b - 4*c - 2*b*alpha - 3*pow(alpha, 2);
        const T aux   = sqrt(-discr);

        x1 = {alpha, 0.0};
        x2 = {(-b - alpha) * 0.5, -aux * 0.5};
        x3 = {(-b - alpha) * 0.5,  aux * 0.5};
    }
    else
    {
        const T pi    = 3.14159265359;
        const T delta = sqrt((b*b - 3*c)/9.0);
        const T h     = 2*delta*delta*delta;
        const T theta = acos(-yn/h)/3;

        const T alpha = xn + 2*delta*cos(theta);
        const T beta  = xn + 2*delta*cos(2*pi/3 - theta);
        const T gamma = xn + 2*delta*cos(2*pi/3 + theta);

        x1 = {alpha, 0.0};
        x2 = {beta,  0.0};
        x3 = {gamma, 0.0};
    }

    return {x1, x2, x3};
}

/// Calculate the root of a non-linear function using Newton's method.
/// @param f The function that returns a pair of \f$ f(x) \f$ and \f$ f^{\prime}(x) \f$.
/// @param x0 The initial guess for the iterative root calculation.
/// @param epsilon The tolerance used in \f$ |f(x)| < \epsilon \f$ to check convergence.
/// @param maxiter The maximum number of iterations.
/// @return A tuple containing the root, the number of iterations, and a boolean flag indicating success if true.
template<typename T>
auto newton(const std::function<std::tuple<T, T>(const T&)>& f, const T& x0, const T& epsilon, std::size_t maxiter) -> std::tuple<T, std::size_t, bool>
{
    using std::abs;
    assert(epsilon > 0.0);
    assert(maxiter > 0);
    T x = x0;
    for(auto i = 0; i < maxiter; ++i)
    {
        const auto [fx, dfx] = f(x);
        x -= fx/dfx;
        if(abs(fx) < epsilon)
            return { x, i, true };
    }
    return { x, maxiter, false };
}

/// Return all real roots of a group of roots
/// @param roots CubicRoots with of complex and real roots
/// @return A vector with all real roots
template<typename T>
auto realRoots(const CubicRoots<T>& roots) -> std::vector<T>
{
    std::vector<T> real_roots;
    if(std::get<0>(roots).imag() == 0.0)
        real_roots.push_back(std::get<0>(roots).real());
    if(std::get<1>(roots).imag() == 0.0)
        real_roots.push_back(std::get<1>(roots).real());
    if(std::get<2>(roots).imag() == 0.0)
        real_roots.push_back(std::get<2>(roots).real());

    return real_roots;
}

} // namespace Reaktoro
