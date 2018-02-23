// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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

#include "Roots.hpp"

// Eigen includes
#include <Reaktoro/Math/Eigen/LU>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>

namespace Reaktoro {

auto cardano(double a, double b, double c, double d) -> CubicRoots
{
    const double xn = -b/(3*a);
    const double yn = a*xn*xn*xn + b*xn*xn + c*xn + d;

    const double d2    = (b*b - 3*a*c)/(9*a*a);
    const double h2    = 4*a*a*d2*d2*d2;
    const double Delta = yn*yn - h2;

    std::complex<double> x1, x2, x3;

    if(Delta > 0.0)
    {
        const double sqrtDelta = std::sqrt(Delta);

        const double operand1 = (-yn + sqrtDelta)/(2*a);
        const double operand2 = (-yn - sqrtDelta)/(2*a);

        const double alpha = xn +
            std::pow(std::abs(operand1), 1.0/3) * std::abs(operand1)/operand1 +
            std::pow(std::abs(operand2), 1.0/3) * std::abs(operand2)/operand2;

        const double discr = b*b - 4*a*c - 2*a*b*alpha - 3*std::pow(a*alpha, 2);
        const double aux   = std::sqrt(-discr);

        x1 = {alpha, 0.0};
        x2 = {(-b - a*alpha)/(2*a), -aux/(2*a)};
        x3 = {(-b - a*alpha)/(2*a),  aux/(2*a)};
    }
    else
    {
        const double pi    = 3.14159265359;
        const double delta = std::sqrt((b*b - 3*a*c)/(9*a*a));
        const double h     = 2*a*delta*delta*delta;
        const double theta = acos(-yn/h)/3;

        const double alpha = xn + 2*delta*cos(theta);
        const double beta  = xn + 2*delta*cos(2*pi/3 - theta);
        const double gamma = xn + 2*delta*cos(2*pi/3 + theta);

        x1 = {alpha, 0.0};
        x2 = {beta,  0.0};
        x3 = {gamma, 0.0};
    }

    return std::make_tuple(x1, x2, x3);
}

auto newton(const std::function<std::tuple<double,double>(double)>& f,
            double x0, double epsilon, unsigned maxiter) -> double
{
    Assert(epsilon > 0.0, "Could not start Newton's method with given parameter.",
        "Expecting a positive tolerance parameter.");
    Assert(maxiter > 0, "Could not start Newton's method with given parameter.",
        "Expecting a positive maximum number of iterations.");
    double x = x0;
    for(unsigned i = 0; i < maxiter; ++i)
    {
        double fx, dfx;
        std::tie(fx, dfx) = f(x);
        x -= fx/dfx;
        if(std::abs(fx) < epsilon)
            return x;
    }
    RuntimeError("Could not find the root of the given non-linear function.",
        "The maximum number of iterations " + std::to_string(maxiter) + " was achieved.");
    return x;
}

auto newton(const std::function<void(VectorConstRef, VectorRef, MatrixRef)>& f,
    VectorConstRef x0, double epsilon, unsigned maxiter) -> Vector
{
    Assert(epsilon > 0.0, "Could not start Newton's method with given parameter.",
        "Expecting a positive tolerance parameter.");
    Assert(maxiter > 0, "Could not start Newton's method with given parameter.",
        "Expecting a positive maximum number of iterations.");
    Index n = x0.rows();
    Vector x = x0;
    Vector fx(n);
    Matrix Jx(n, n);
    for(unsigned i = 0; i < maxiter; ++i)
    {
        f(x, fx, Jx);
        x -= Jx.lu().solve(fx);
        if(max(abs(fx)) < epsilon)
            return x;
    }
    RuntimeError("Could not find the root of the given non-linear function.",
        "The maximum number of iterations " + std::to_string(maxiter) + " was achieved.");
    return x;
}

} // namespace Reaktoro
