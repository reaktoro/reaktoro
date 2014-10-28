// Reaktor is a C++ library for computational reaction modelling.
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

#include "AlgorithmUtils.hpp"

// C++ includes
#include <limits>

namespace Reaktor {

auto largestStep(const Vector& p, const Vector& dp) -> double
{
    Vector res = -p/dp;
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

} // namespace Reaktor
