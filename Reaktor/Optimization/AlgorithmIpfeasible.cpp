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
    const unsigned n = problem.numVariables();
    const unsigned m = problem.numConstraints();
    const unsigned t = n + 2*m;

    const double mu = options.ipnewton.mu;
    const double rho = std::sqrt(options.ipnewton.mu);

    ObjectiveResult objres;
    objres.grad.resize(t);
    objres.grad << zeros(n), ones(2*m);

    objres.hessian = zeros(t, t);
    objres.hessian.block(0, 0, n, n) = rho * identity(n, n);

    const Vector xr = (result.solution.x.size() == n) ? result.solution.x : mu * ones(n);

    result.solution.x = xr;

    ObjectiveFunction objective = [=](const Vector& x) mutable -> ObjectiveResult
    {
        const Vector xx = rows(x, 0, n);
        const Vector xp = rows(x, n, m);
        const Vector xn = rows(x, n + m, m);
        objres.func = (xp + xn).sum() + 0.5 * rho * (xx - xr).dot(xx - xr);
        rows(objres.grad, 0, n) = rho*(xx - xr);
        return objres;
    };

    ConstraintFunction constraint = [=](const Vector& x) mutable -> ConstraintResult
    {
        auto xx = rows(x, 0, n);
        auto xp = rows(x, n, m);
        auto xn = rows(x, n + m, m);

        ConstraintResult res = problem.constraint()(xx);
        res.func += xn - xp;
        res.grad = arma::join_horiz(res.grad, -arma::eye(m, m));
        res.grad = arma::join_horiz(res.grad,  arma::eye(m, m));

        return res;
    };

    options.ipnewton.saddle_point.algorithm = Rangespace;
    options.ipnewton.saddle_point.properties = DiagonalH;

    problem = OptimumProblem(t, m);
    problem.setObjective(objective);
    problem.setConstraint(constraint);
    problem.setLowerBounds(0);

    Vector xx = result.solution.x;
    Vector xp = arma::zeros(m) + options.ipnewton.mu;
    Vector xn = arma::zeros(m) + options.ipnewton.mu;

    result.solution.x  = arma::join_vert(arma::join_vert(xx, xp), xn);
    result.solution.y  = arma::zeros(m);
    result.solution.zl = options.ipnewton.mu/result.solution.x;

    ipnewton(problem, result, options);

    result.solution.x  = result.solution.x.rows(0, n-1);
    result.solution.y  = arma::zeros(m);
    result.solution.zl = result.solution.x/options.ipnewton.mu;
}

} // namespace Reaktor
