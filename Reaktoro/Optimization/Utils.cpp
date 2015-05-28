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

#include "Utils.hpp"

// C++ includes
#include <limits>

namespace Reaktoro {

auto largestStep(const Vector& p, const Vector& dp) -> double
{
    Vector res = -p.array() / dp.array();
    double alpha = infinity();
    for(unsigned i = 0; i < res.size(); ++i)
        if(res[i] > 0.0 and res[i] < alpha)
            alpha = res[i];
    return alpha;
}

auto fractionToTheBoundary(const Vector& p, const Vector& dp, double tau) -> double
{
    double alpha_max = 1.0;
    for(unsigned i = 0; i < p.size(); ++i)
        if(dp[i] < 0.0) alpha_max = std::min(alpha_max, -tau*p[i]/dp[i]);
    return alpha_max;
}

auto lessThan(double lhs, double rhs, double baseval) -> bool
{
    const double epsilon = std::numeric_limits<double>::epsilon();
    return lhs < rhs + 10.0 * epsilon * std::abs(baseval);
}

auto greaterThan(double lhs, double rhs, double baseval) -> bool
{
    const double epsilon = std::numeric_limits<double>::epsilon();
    return lhs > rhs - 10.0 * epsilon * std::abs(baseval);
}

auto infinity() -> double
{
    return std::numeric_limits<double>::infinity();
}

auto bfgs() -> std::function<Matrix(const Vector&, const Vector&)>
{
    Vector x0;
    Vector g0;
    Vector dx;
    Vector dg;
    Matrix H;

    std::function<Matrix(const Vector&, const Vector&)> f = [=](const Vector& x, const Vector& g) mutable
    {
        if(x0.size() == 0)
        {
            x0.noalias() = x;
            g0.noalias() = g;
            H = diag(x);
            return H;
        }

        dx.noalias() = x - x0;
        dg.noalias() = g - g0;
        x0.noalias() = x;
        g0.noalias() = g;

        const unsigned n = x.size();
        const double a = dx.dot(dg);
        const auto I = identity(n, n);

        H = (I - dx*tr(dg)/a)*H*(I - dg*tr(dx)/a) + dx*tr(dx)/a;

        return H;
    };

    return f;
}

auto minimizeGoldenSectionSearch(const std::function<double(double)>& f, double a, double b, double tol) -> double
{
    //---------------------------------------------------------------
    // Reference: http://en.wikipedia.org/wiki/Golden_section_search
    //---------------------------------------------------------------

    // The golden ratio
    const double phi = 0.61803398875;

    double c = b - phi*(b - a);
    double d = a + phi*(b - a);

    while(std::abs(c - d) > tol*(std::abs(a) + std::abs(b)))
    {
        double fc = f(c);
        double fd = f(d);

        if(fc < fd)
        {
            b = d;
            d = c; //  #fd=fc;fc=f(c)
            c = b - phi*(b - a);
        }
        else
        {
            a = c;
            c = d; // fc=fd;fd=f(d)
            d = a + phi*(b - a);
        }
    }

    return (b + a)/2.0;
}

} // namespace Reaktoro
