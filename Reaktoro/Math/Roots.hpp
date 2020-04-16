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
    using std::cbrt;
    using std::cos;
    using std::sqrt;
    using std::complex;

    const auto xn = -b/3;
    const auto yn = d + xn*(c + xn*(b + xn));

    const auto delta2 = (b*b - 3*c)/9;

    const auto u = yn*yn;
    const auto v = 4*delta2*delta2*delta2;

    const auto discr = u - v;

    const auto eps = 100*std::numeric_limits<T>::epsilon();

    if(discr > eps) {
        const auto sqrtdiscr = sqrt(discr);
        const auto aux1 = cbrt( 0.5*(-yn + sqrtdiscr) );
        const auto aux2 = cbrt( 0.5*(-yn - sqrtdiscr) );
        const auto alpha = xn + aux1 + aux2;
        const auto z1 = alpha - xn;
        const auto aux3 = 0.5*sqrt(3*(z1*z1 - 4*delta2));
        const auto aux4 = xn - 0.5*z1;
        const auto beta = complex{aux4, -aux3};
        const auto gamma = complex{aux4, aux3};
        return {alpha, beta, gamma};
    }
    else if(discr < -eps) {
        const auto pi = 3.14159265358979323846;
        const auto delta = sqrt(delta2);
        const auto h = 2*delta2*delta;
        const auto theta = acos(-yn/h) / 3.0;
        const auto alpha = xn + 2*delta*cos(theta);
        const auto beta  = xn + 2*delta*cos(2*pi/3 - theta);
        const auto gamma = xn + 2*delta*cos(2*pi/3 + theta);
        return {gamma, beta, alpha};
    }
    else {
        const auto delta = cbrt( 0.5*yn );
        const auto alpha = xn + delta;
        const auto gamma = xn - 2*delta;
        if(alpha < gamma) return {alpha, alpha, gamma};
        else return {gamma, alpha, alpha};
    }
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
        assert(dfx != 0.0);
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
