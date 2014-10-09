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
#include <algorithm>

// Reaktor includes
#include <Reaktor/Common/Macros.hpp>

namespace Reaktor {

OptimumProblem::OptimumProblem(unsigned n, unsigned m)
: n(n), m(m), l(arma::zeros(n)), u(INFINITY*arma::ones(n))
{}

auto OptimumProblem::setObjective(const ObjectiveFunction& objective) -> void
{
    f = objective;
}

auto OptimumProblem::setConstraint(const ConstraintFunction& constraint) -> void
{
    h = constraint;
}

auto OptimumProblem::setLowerBounds(const Vector& lower) -> void
{
    Assert(lower.size() == n, "Dimension of the upper bound vector does not match dimension of primal variables.");
    l = lower;
}

auto OptimumProblem::setUpperBounds(const Vector& upper) -> void
{
    Assert(upper.size() == n, "Dimension of the upper bound vector does not match dimension of primal variables.");
    u = upper;
}

auto OptimumProblem::numVariables() const -> unsigned
{
    return n;
}

auto OptimumProblem::numConstraints() const -> unsigned
{
    return m;
}

auto OptimumProblem::objective() const -> const ObjectiveFunction&
{
    return f;
}

auto OptimumProblem::constraint() const -> const ConstraintFunction&
{
    return h;
}

auto OptimumProblem::lowerBounds() const -> const Vector&
{
    return l;
}

auto OptimumProblem::upperBounds() const -> const Vector&
{
    return u;
}

auto dominated(const FilterEntry& a, const FilterEntry& b) -> bool
{
    for(unsigned i = 0; i < a.size(); ++i)
        if(a[i] < b[i]) return false;
    return true;
}

auto acceptable(const Filter& filter, const FilterEntry& entry) -> bool
{
    for(const FilterEntry& point : filter)
        if(dominated(entry, point)) return false;
    return true;
}

auto extend(Filter& filter, const FilterEntry& entry) -> void
{
    // Check if the entry trying to be added to the filter is accepted to it
    if(acceptable(filter, entry))
    {
        // Define the domination function to remove dominated points from the filter
        auto is_dominated = [=](const FilterEntry& point)
        {
            return dominated(point, entry);
        };

        // Remove all dominated entries in the filter
        filter.erase(std::remove_if(filter.begin(), filter.end(), is_dominated), filter.end());

        // Add the new entry to the filter
        filter.push_back(entry);
    }
}

auto largestStep(const Vector& p, const Vector& dp) -> double
{
    Vector res = -p/dp;
    double alpha = INFINITY;
    for(unsigned i = 0; i < res.size(); ++i)
        if(res[i] > 0.0 and res[i] < alpha)
            alpha = res[i];
    return alpha;
}

auto fractionToTheBoundary(const Vector& p, const Vector& dp, double tau) -> double
{
    double alpha_max = 1.0;
    for(unsigned i = 0; i < p.size(); ++i)
    {
        const double alpha = -tau*p[i]/dp[i];
        if(0.0 < alpha and alpha < alpha_max)
            alpha_max = alpha;
    }
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

} // namespace Reaktor
