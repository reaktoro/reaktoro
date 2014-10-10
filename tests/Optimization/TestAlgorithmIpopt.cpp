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

#include "TestAlgorithmIpopt.hpp"

// Reaktor includes
#include <Reaktor/Reaktor.hpp>

namespace Reaktor {
namespace {

auto test_ipfeasible() -> void
{
    const unsigned n = 5;
    const unsigned m = 2;

    OptimumProblem problem(n, m);
    OptimumResult result;
    OptimumOptions options;

    ipfeasible(problem, result, options);

    ASSERT_EQUAL(n, result.solution.x.size());
    ASSERT_EQUAL(m, result.solution.y.size());
    ASSERT_EQUAL(n, result.solution.zl.size());
    ASSERT_EQUAL(n, result.solution.zu.size());
    ASSERT(arma::all(result.solution.x > problem.lowerBounds()));
    ASSERT(arma::all(result.solution.zl == options.mu));
    ASSERT(arma::all(result.solution.zu == 0.0));
}


} // namespace


auto testSuiteAlgorithmIpopt() -> cute::suite
{
    cute::suite s;

    s += CUTE(test_ipfeasible);

    return s;
}

} // namespace Reaktor
