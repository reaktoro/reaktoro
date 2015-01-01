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

#include "OptimumSolverIpfeasible.hpp"

// Reaktor includes
#include <Reaktor/Optimization/OptimumOptions.hpp>
#include <Reaktor/Optimization/OptimumProblem.hpp>
#include <Reaktor/Optimization/OptimumResult.hpp>
#include <Reaktor/Optimization/OptimumSolverIpnewton.hpp>

namespace Reaktor {

struct OptimumSolverIpfeasible::Impl
{
    OptimumSolverIpnewton ipnewton;

    auto approximate(OptimumProblem problem, OptimumResult& result, OptimumOptions options) -> void;
};

auto OptimumSolverIpfeasible::Impl::approximate(OptimumProblem problem, OptimumResult& result, OptimumOptions options) -> void
{
    // Auxiliary variables
    const unsigned n = problem.numVariables();
    const unsigned m = problem.numConstraints();
    const double  mu = options.ipnewton.mu;
    const double rho = std::sqrt(mu);

    // The total number of primal variables for the feasibility problem
    const unsigned t = n + 2*m;

    // Initialize the solution variables if it has not been done before
    if(result.solution.x.size() != n)
        result.solution.x = mu * ones(n);

    // The reference point from which the solution cannot differ much
    const Vector xr = result.solution.x;

    // Define the objective function of the feasibility problem
    ObjectiveFunction objective = [=](const Vector& x) mutable
    {
        const auto xx = rows(x, 0, n);
        const auto xp = rows(x, n, m);
        const auto xn = rows(x, n + m, m);
        return (xp + xn).sum() + 0.5 * rho * (xx - xr).dot(xx - xr);
    };

    // Define the gradient function of the objective function of the feasibility problem
    Vector g(t);
    g << zeros(n), ones(2*m);
    ObjectiveGradFunction objective_grad = [=](const Vector& x) mutable
    {
        const auto xx = rows(x, 0, n);
        rows(g, 0, n) = rho*(xx - xr);
        return g;
    };

    // Define the Hessian function of the objective function of the feasibility problem
    Vector diagH = zeros(t);
    rows(diagH, 0, n) = rho * ones(n);
    ObjectiveDiagonalHessianFunction objective_hessian = [=](const Vector& x) mutable
    {
        return diagH;
    };

    // Define the equality constraint function of the feasibility problem
    Vector h;
    ConstraintFunction constraint = [=](const Vector& x) mutable
    {
        const auto xx = rows(x, 0, n);
        const auto xp = rows(x, n, m);
        const auto xn = rows(x, n + m, m);
        h  = problem.constraint(xx);
        h += xn - xp;
        return h;
    };

    // Define the gradient function of the equality constraint function of the feasibility problem
    Vector A(m, t);
    cols(A, n, m)     = -identity(m, m);
    cols(A, n + m, m) =  identity(m, m);
    ConstraintGradFunction constraint_grad = [=](const Vector& x) mutable
    {
        const auto xx = rows(x, 0, n);
        cols(A, 0, n) = problem.constraintGrad(xx);
        return A;
    };

    // Set the initial guess
    const Vector xx = result.solution.x;
    const Vector xp = mu * ones(m);
    const Vector xn = mu * ones(m);

    result.solution.x.resize(t);
    result.solution.x << xx, xp, xn;
    result.solution.y  = zeros(m);
    result.solution.z = mu/result.solution.x.array();

    // Define the feasibility problem
    problem = OptimumProblem(t, m);
    problem.setObjective(objective);
    problem.setObjectiveGrad(objective_grad);
    problem.setObjectiveHessian(objective_hessian);
    problem.setConstraint(constraint);
    problem.setLowerBounds(0);

    // Solve the feasibility problem
    ipnewton.solve(problem, result, options);

    // Prepare the exported result of the calculation
    result.solution.x  = rows(result.solution.x, 0, n);
    result.solution.y  = zeros(m);
    result.solution.z = mu/result.solution.x.array();
}

OptimumSolverIpfeasible::OptimumSolverIpfeasible()
: pimpl(new Impl())
{}

OptimumSolverIpfeasible::OptimumSolverIpfeasible(const OptimumSolverIpfeasible& other)
: pimpl(new Impl(*other.pimpl))
{}

OptimumSolverIpfeasible::~OptimumSolverIpfeasible()
{}

auto OptimumSolverIpfeasible::operator=(OptimumSolverIpfeasible other) -> OptimumSolverIpfeasible&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto OptimumSolverIpfeasible::approximate(const OptimumProblem& problem, OptimumResult& result) -> void
{
    approximate(problem, result, {});
}

auto OptimumSolverIpfeasible::approximate(const OptimumProblem& problem, OptimumResult& result, const OptimumOptions& options) -> void
{
    pimpl->approximate(problem, result, options);
}

} // namespace Reaktor
