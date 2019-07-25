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

#include "Roots.hpp"

// Eigen includes
#include <Reaktoro/deps/eigen3/Eigen/LU>

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

auto ferrari(const double a, const double b, const double c, const double d, const double e)->QuarticRoots
{
    const double Delta = 256.0 * std::pow(a*e, 3.0) - 192.0 * std::pow(a*e, 2.0) * b * d
                        - 128.0 * std::pow(a*c*e, 2.0) + 144.0 * std::pow(a*d, 2.0) * c*e 
                        - 27.0 * std::pow(a*std::pow(d, 2.0), 2.0) + 144.0 * std::pow(b*e, 2.0) * a*c 
                        - 6.0 * std::pow(b*d, 2.0) * a*e - 80.0 * a*b*std::pow(c, 2.0) * d*e 
                        + 18.0 * a*b*c*std::pow(d, 3.0) + 16.0 * a*std::pow(c, 4.0) * e 
                        - 4.0 * a*std::pow(c, 3.0) * std::pow(d, 2.0) - 27.0 * std::pow(b, 4.0) * std::pow(e, 2.0) 
                        + 18.0 * std::pow(b, 3.0) * c*d*e - 4.0 * std::pow(b*d, 3.0) 
                        - 4.0 * std::pow(b, 2.0) * std::pow(c, 3.0) * e + std::pow(b*c*d, 2.0);

    const double Delta0 = std::pow(c, 2.0) - 3.0 * b*d + 12.0 * a*e;
    const double P = 8.0 * a*c - 3.0 * std::pow(b, 2.0);
    const double R = std::pow(b, 3.0) + 8.0 * d*std::pow(a, 2.0) - 4 * a*b*c;
    const double D = 64.0 * std::pow(a, 3.0) * e - 16.0 * std::pow(a*c, 2.0) 
                    + 16.0 * a*std::pow(b, 2.0) * c - 16 * std::pow(a, 2.0) * b*d - 3 * std::pow(b, 4.0);
    
    //reduced to a cubic equation
    const double A = c / a - (3.0 / 8.0)*std::pow(b / a, 2.0);
    const double B = d / a - b * c / (2.0 * std::pow(a, 2.0)) + 1.0 / 8.0 * std::pow(b / a, 3.0);
    const double C = e / a - b * d / (4.0 * std::pow(a, 2.0)) + std::pow(b, 2.0) * c / (16.0 * std::pow(a, 3.0)) - 3.0 / 256.0 * std::pow(b / a, 4.0);

    const double coef1 = 8.0;
    const double coef2 = -4.0 * A;
    const double coef3 = -8.0 * C;
    const double coef4 = 4.0 * A * C - std::pow(B, 2.0);

    auto cubic_roots = cardano(coef1, coef2, coef3, coef4);
    auto cubic_r1 = std::get<0>(cubic_roots);

    //compute quartic roots from cubic roots
    auto y1 = -(1.0 / 2.0)*std::pow(2.0 * cubic_r1 - A, 0.5) + (1.0 / 2.0)*std::pow(-2.0 * cubic_r1 - A + 2.0 * B / (std::sqrt(2.0 * cubic_r1 - A)), 0.5);
    auto y2 = -(1.0 / 2.0)*std::pow(2.0 * cubic_r1 - A, 0.5) - (1.0 / 2.0)*std::pow(-2.0 * cubic_r1 - A + 2.0 * B / (std::sqrt(2.0 * cubic_r1 - A)), 0.5);
    auto y3 = (1.0 / 2.0)*std::pow(2.0 * cubic_r1 - A, 0.5) + (1.0 / 2.0)*std::pow(-2.0 * cubic_r1 - A - 2.0 * B / (std::sqrt(2.0 * cubic_r1 - A)), 0.5);
    auto y4 = (1.0 / 2.0)*std::pow(2.0 * cubic_r1 - A, 0.5) - (1.0 / 2.0)*std::pow(-2.0 * cubic_r1 - A - 2.0 * B / (std::sqrt(2.0 * cubic_r1 - A)), 0.5);

    auto x1 = y1 - b / (4.0 * a);
    auto x2 = y2 - b / (4.0 * a);
    auto x3 = y3 - b / (4.0 * a);
    auto x4 = y4 - b / (4.0 * a);

    return std::make_tuple(x1, x2, x3, x4);
}

auto realRoots(const QuarticRoots& roots) -> std::vector<double>
{
    std::vector<double> real_roots;
    if (std::get<0>(roots).imag() == 0)
        real_roots.push_back(std::get<0>(roots).real());
    if (std::get<1>(roots).imag() == 0)
        real_roots.push_back(std::get<1>(roots).real());
    if (std::get<2>(roots).imag() == 0)
        real_roots.push_back(std::get<2>(roots).real());
    if (std::get<3>(roots).imag() == 0)
        real_roots.push_back(std::get<3>(roots).real());
    
    return real_roots;
}

auto realRoots(const CubicRoots& roots) -> std::vector<double>
{
    std::vector<double> real_roots;
    if (std::get<0>(roots).imag() == 0)
        real_roots.push_back(std::get<0>(roots).real());
    if (std::get<1>(roots).imag() == 0)
        real_roots.push_back(std::get<1>(roots).real());
    if (std::get<2>(roots).imag() == 0)
        real_roots.push_back(std::get<2>(roots).real());

    return real_roots;
}

} // namespace Reaktoro
