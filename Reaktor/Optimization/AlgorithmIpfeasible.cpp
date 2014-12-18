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

#include "AlgorithmIpfeasible.hpp"

// Reaktor includes
#include <Reaktor/Optimization/OptimumOptions.hpp>
#include <Reaktor/Optimization/OptimumProblem.hpp>
#include <Reaktor/Optimization/OptimumResult.hpp>
#include <Reaktor/Optimization/AlgorithmIpnewton.hpp>

namespace Reaktor {

auto ipfeasible(OptimumProblem problem, OptimumResult& result, OptimumOptions options) -> void
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

	// The result of the objective function evaluation
    ObjectiveResult objective_result;
    objective_result.grad.resize(t);
    objective_result.grad << zeros(n), ones(2*m);
    objective_result.hessian = zeros(t, t);
    objective_result.hessian.block(0, 0, n, n) = rho * identity(n, n);

    // Define the objective function of the feasibility problem
    ObjectiveFunction objective = [=](const Vector& x) mutable -> ObjectiveResult
    {
        const auto xx = rows(x, 0, n);
        const auto xp = rows(x, n, m);
        const auto xn = rows(x, n + m, m);
        objective_result.func = (xp + xn).sum() + 0.5 * rho * (xx - xr).dot(xx - xr);
        rows(objective_result.grad, 0, n) = rho*(xx - xr);
        return objective_result;
    };

	// The result of the constraint function evaluation
    ConstraintResult constraint_result, temp;
    constraint_result.grad.resize(m, t);
    cols(constraint_result.grad, n    , m) = -identity(m, m);
    cols(constraint_result.grad, n + m, m) =  identity(m, m);

    // Define the constraint function of the feasibility problem
    ConstraintFunction constraint = [=](const Vector& x) mutable -> ConstraintResult
    {
        const auto xx = rows(x, 0, n);
        const auto xp = rows(x, n, m);
        const auto xn = rows(x, n + m, m);
        temp = problem.constraint()(xx);
        constraint_result.func = temp.func + xn - xp;
        cols(constraint_result.grad, 0, n) = temp.grad;
        return constraint_result;
    };

    // Set some options for efficient calculation of the saddle point problems (KKT equations)
    options.ipnewton.saddle_point.algorithm = Rangespace;
    options.ipnewton.saddle_point.properties = DiagonalH;

    // Set the initial guess
    const Vector xx = result.solution.x;
    const Vector xp = mu * ones(m);
    const Vector xn = mu * ones(m);

    result.solution.x.resize(t);
    result.solution.x << xx, xp, xn;
    result.solution.y  = zeros(m);
    result.solution.zl = mu/result.solution.x.array();

    // Define the feasibility problem
	problem = OptimumProblem(t, m);
    problem.setObjective(objective);
    problem.setConstraint(constraint);
    problem.setLowerBounds(0);

    // Solve the feasibility problem
    ipnewton(problem, result, options);

    // Prepare the exported result of the calculation
    result.solution.x  = rows(result.solution.x, 0, n);
    result.solution.y  = zeros(m);
    result.solution.zl = mu/result.solution.x.array();
}

} // namespace Reaktor
