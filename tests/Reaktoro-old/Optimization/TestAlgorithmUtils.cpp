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

#include "TestAlgorithmUtils.hpp"

// Reaktoro includes
#include <Reaktoro/Reaktoro.hpp>

namespace Reaktoro {
namespace {

#define ASSERT_EQUAL_ARMA(expected, actual) ASSERT(arma::all(expected==actual))

auto test_OptimumProblem() -> void
{
    const unsigned m = 2;
    const unsigned n = 5;

    ObjectiveFunction objective = [](const auto& x)
    {
        ObjectiveResult f;
        f.func = arma::sum(x);
        f.grad = arma::ones(n);
        f.hessian = arma::zeros(n, n);
        return f;
    };

    ConstraintFunction constraint = [](const auto& x)
    {
        Matrix A = arma::zeros(m, n);
        A.submat(0, 0, m-1, m-1) = arma::eye(m, m);
        Vector b = arma::ones(m);
        ConstraintResult h;
        h.func = A*x-b;
        h.grad = A;
        h.hessian = arma::zeros(m, n, n);
        return h;
    };

    OptimumProblem problem(n, m);
    problem.setObjective(objective);
    problem.setConstraint(constraint);

    ASSERT_EQUAL(n, problem.numVariables());
    ASSERT_EQUAL(m, problem.numConstraints());
    ASSERT_EQUAL_ARMA(arma::zeros(n), problem.lowerBounds());
    ASSERT_EQUAL_ARMA(INFINITY*arma::ones(n), problem.upperBounds());

    problem.setLowerBounds(arma::ones(n) * 2.0);
    problem.setUpperBounds(arma::ones(n) * 5.0);

    ASSERT_EQUAL_ARMA(arma::ones(n) * 2.0, problem.lowerBounds());
    ASSERT_EQUAL_ARMA(arma::ones(n) * 5.0, problem.upperBounds());
}

auto test_dominated() -> void
{
    ASSERT(dominated({2.0, 2.0}, {1.0, 1.0}));
    ASSERT(dominated({1.0, 3.0}, {0.9, 2.0}));
    ASSERT(dominated({2.0, 3.0}, {2.0, 3.0}));
    ASSERT(not dominated({1.0, 1.0}, {2.0, 3.0}));
}

auto test_acceptable() -> void
{
    Filter filter;
    filter.push_back({1.0, 4.0});
    filter.push_back({2.0, 3.0});
    filter.push_back({3.0, 2.0});
    filter.push_back({4.0, 1.0});

    ASSERT(acceptable(filter, {0.0, 0.0}));
    ASSERT(acceptable(filter, {1.0, 1.0}));
    ASSERT(acceptable(filter, {2.0, 2.0}));
    ASSERT(acceptable(filter, {1.0, 2.0}));
    ASSERT(acceptable(filter, {2.0, 1.0}));
    ASSERT(acceptable(filter, {2.0, 2.0}));
    ASSERT(acceptable(filter, {1.0, 3.0}));
    ASSERT(not acceptable(filter, {1.0, 4.0}));
    ASSERT(not acceptable(filter, {1.0, 5.0}));
}

auto test_extend() -> void
{
    Filter filter;
    filter.push_back({1.0, 4.0});
    filter.push_back({2.0, 3.0});
    filter.push_back({3.0, 2.0});
    filter.push_back({4.0, 1.0});

    Filter filter1 = filter;
    extend(filter1, {2.0, 2.0});
    Filter expected1 = {{1.0, 4.0}, {2.0, 2.0}, {4.0, 1.0}};

    Filter filter2 = filter;
    extend(filter2, {5.0, 5.0});
    Filter expected2 = filter;

    Filter filter3 = filter;
    extend(filter3, {2.5, 0.5});
    Filter expected3 = {{1.0, 4.0}, {2.0, 3.0}, {2.5, 0.5}};

    ASSERT(equal(expected1, filter1));
    ASSERT(equal(expected2, filter2));
    ASSERT(equal(expected3, filter3));
}

auto test_largestStep() -> void
{
    Vector p   = {1.0, 1.0};
    Vector dp1 = {1.0, 1.0};
    Vector dp2 = {-2.0, 1.0};
    Vector dp3 = {1.0, -5.0};
    Vector dp4 = {-1.0, -0.5};

    ASSERT(largestStep(p, dp1) == INFINITY);
    ASSERT_EQUAL(0.5, largestStep(p, dp2));
    ASSERT_EQUAL(0.2, largestStep(p, dp3));
    ASSERT_EQUAL(1.0, largestStep(p, dp4));
}

auto test_fractionToTheBoundary() -> void
{
    Vector p   = {1.0, 1.0};
    Vector dp1 = {1.0, 1.0};
    Vector dp2 = {-2.0, 1.0};
    Vector dp3 = {1.0, -5.0};
    Vector dp4 = {-1.0, -0.5};

    const double tau = 0.995;

    ASSERT_EQUAL(1.0, fractionToTheBoundary(p, dp1, tau));
    ASSERT_EQUAL(tau * 0.5, fractionToTheBoundary(p, dp2, tau));
    ASSERT_EQUAL(tau * 0.2, fractionToTheBoundary(p, dp3, tau));
    ASSERT_EQUAL(tau * 1.0, fractionToTheBoundary(p, dp4, tau));
}

auto test_lessThan() -> void
{
    ASSERT(lessThan(1.0, 2.0, 1.0));
    ASSERT(lessThan(1e-6, 1e-6-1e-15, 1.0));
    ASSERT(not lessThan(1e-6+2e-15, 1e-6-1e-15, 1.0));
    ASSERT(not lessThan(2.0, 1.0, 1.0));
}

auto test_greaterThan() -> void
{
    ASSERT(not greaterThan(1.0, 2.0, 1.0));
    ASSERT(greaterThan(1e-6, 1e-6+1e-15, 1.0));
    ASSERT(greaterThan(1e-6+2e-15, 1e-6-1e-15, 1.0));
    ASSERT(greaterThan(2.0, 1.0, 1.0));
}

} // namespace

auto testSuiteAlgorithmUtils() -> cute::suite
{
    cute::suite s;

    s += CUTE(test_OptimumProblem);
    s += CUTE(test_dominated);
    s += CUTE(test_acceptable);
    s += CUTE(test_extend);
    s += CUTE(test_largestStep);
    s += CUTE(test_fractionToTheBoundary);
    s += CUTE(test_lessThan);
    s += CUTE(test_greaterThan);

    return s;
}

} // namespace Reaktoro
