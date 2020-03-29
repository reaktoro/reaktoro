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

#include "Utils.hpp"

// C++ includes
#include <limits>

namespace Reaktoro {

auto largestStep(VectorXdConstRef p, VectorXdConstRef dp) -> double
{
    VectorXd res = -p.array() / dp.array();
    double alpha = infinity();
    for(unsigned i = 0; i < res.size(); ++i)
        if(res[i] > 0.0 && res[i] < alpha)
            alpha = res[i];
    return alpha;
}

auto fractionToTheBoundary(VectorXdConstRef p, VectorXdConstRef dp, double tau) -> double
{
    Index i;
    return fractionToTheBoundary(p, dp, tau, i);
}

auto fractionToTheBoundary(VectorXdConstRef p, VectorXdConstRef dp, double tau, Index& ilimiting) -> double
{
    ilimiting = p.size();
    double alpha_max = 1.0;
    for(int i = 0; i < p.size(); ++i)
    {
        if(dp[i] < 0.0)
        {
            const double alpha_trial = -tau*p[i]/dp[i];
            if(alpha_trial < alpha_max)
            {
                alpha_max = alpha_trial;
                ilimiting = i;
            }
        }
    }

    return alpha_max;
}

auto fractionToTheBoundary(VectorXdConstRef p, VectorXdConstRef dp, MatrixXdConstRef C, VectorXdConstRef r, double tau) -> double
{
    // The number of linear inequality constraints
    const Index m = C.rows();

    // Check if there is any inequality constraint and return 1.0 if not
    if(m == 0) return 1.0;

    // Otherwise, compute max(alpha)
    double alpha_max = 1.0;
    for(Index i = 0; i < m; ++i)
    {
        const double tmp = C.row(i).dot(dp);
        if(tmp < 0.0)
        {
            const double alpha_trial = -tau*(C.row(i).dot(p) - r[i])/tmp;
            if(alpha_trial < alpha_max)
                alpha_max = alpha_trial;
        }
    }

    return alpha_max;
}

auto fractionToTheLowerBoundary(VectorXdConstRef p, VectorXdConstRef dp, VectorXdConstRef lower, double tau) -> double
{
    double alpha_max = 1.0;
    for(unsigned i = 0; i < p.size(); ++i)
        if(dp[i] < 0.0) alpha_max = std::min(alpha_max, -tau*(p[i] - lower[i])/dp[i]);
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

auto bfgs() -> std::function<MatrixXd(VectorXdConstRef, VectorXdConstRef)>
{
    VectorXd x0;
    VectorXd g0;
    VectorXd dx;
    VectorXd dg;
    MatrixXd H;

    std::function<MatrixXd(VectorXdConstRef, VectorXdConstRef)> f = [=](VectorXdConstRef x, VectorXdConstRef g) mutable
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

auto minimizeGoldenSectionSearch(const std::function<double(double)>& f, double tol) -> double
{
    //---------------------------------------------------------------
    // Reference: http://en.wikipedia.org/wiki/Golden_section_search
    //---------------------------------------------------------------

    // The golden ratio
    const double phi = 0.61803398875;

    double a = 0.0;
    double b = 1.0;

    double c = 1 - phi;
    double d = phi;

    if(std::abs(c - d) < tol)
        return (b + a)/2.0;

    double fc = f(c);
    double fd = f(d);

    while(std::abs(c - d) > tol)
    {
        if(fc < fd)
        {
            b = d;
            d = c;
            c = b - phi*(b - a);
            fd = fc;
            fc = f(c);
        }
        else
        {
            a = c;
            c = d;
            d = a + phi*(b - a);
            fc = fd;
            fd = f(d);
        }
    }

    return (b + a)/2.0;
}


auto minimizeGoldenSectionSearch(const std::function<double(double)>& f, double a, double b, double tol) -> double
{
    auto g = [=](double x)
    {
        return f(a + x*(b - a));
    };

    const double xmin = minimizeGoldenSectionSearch(g, tol);
    return a + xmin*(b - a);
}

auto minimizeBrent(const std::function<double(double)>& f, double min, double max, double tolerance, unsigned maxiters) -> double
{
    //-------------------------------------------------------------------
    // The code below was adapted from boost library, found at header
    // boost/math/tools/minima.hpp under the name brent_find_minima.
    //-------------------------------------------------------------------
    double x; // minima so far
    double w; // second best point
    double v; // previous value of w
    double u; // most recent evaluation point
    double delta; // The distance moved in the last step
    double delta2; // The distance moved in the step before last
    double fu, fv, fw, fx; // function evaluations at u, v, w, x
    double mid; // midpoint of min and max
    double fract1, fract2; // minimal relative movement in x

    const double golden = 0.3819660;// golden ratio, don't need too much precision here!

    x = w = v = max;
    fw = fv = fx = f(x);
    delta2 = delta = 0;

    unsigned count = maxiters;

    do
    {
        // get midpoint
        mid = (min + max)/2;

        // work out if we're done already:
        fract1 = tolerance * fabs(x) + tolerance/4;
        fract2 = 2 * fract1;

        if(fabs(x - mid) <= (fract2 - (max - min)/2))
            break;

        if(fabs(delta2) > fract1)
        {
            // try and construct a parabolic fit:
            double r = (x - w) * (fx - fv);
            double q = (x - v) * (fx - fw);
            double p = (x - v) * q - (x - w) * r;
            q = 2 * (q - r);
            if(q > 0)
            p = -p;
            q = fabs(q);
            double td = delta2;
            delta2 = delta;

            // determine whether a parabolic step is acceptible or not:
            if((fabs(p) >= fabs(q * td / 2)) || (p <= q * (min - x)) || (p >= q * (max - x)))
            {
                // nope, try golden section instead
                delta2 = (x >= mid) ? min - x : max - x;
                delta = golden * delta2;
            }
            else
            {
                // when, parabolic fit:
                delta = p / q;
                u = x + delta;
                if(((u - min) < fract2) || ((max- u) < fract2))
                delta = (mid - x) < 0 ? -fabs(fract1) : fabs(fract1);
            }
        }
        else
        {
            // golden section:
            delta2 = (x >= mid) ? min - x : max - x;
            delta = golden * delta2;
        }
        // update current position:
        u = (fabs(delta) >= fract1) ? x + delta : (delta > 0 ? x + fabs(fract1) : x - fabs(fract1));
        fu = f(u);
        if(fu <= fx)
        {
            // good new point is an improvement!
            // update brackets:
            if(u >= x)
            min = x;
            else
            max = x;
            // update control points:
            v = w;
            w = x;
            x = u;
            fv = fw;
            fw = fx;
            fx = fu;
        }
        else
        {
            // Oh dear, point u is worse than what we have already,
            // even so it *must* be better than one of our endpoints:
            if(u < x)
            min = u;
            else
            max = u;
            if((fu <= fw) || (w == x))
            {
                // however it is at least second best:
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            }
            else if((fu <= fv) || (v == x) || (v == w))
            {
                // third best:
                v = u;
                fv = fu;
            }
        }

    } while(--count);

    maxiters -= count;

    return x;
}

} // namespace Reaktoro
