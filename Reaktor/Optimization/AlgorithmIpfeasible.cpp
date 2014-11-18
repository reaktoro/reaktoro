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
    objres.hessian = arma::zeros(t, t);
    objres.hessian.submat(0, 0, n-1, n-1) = rho*arma::eye(n, n);
    objres.grad = arma::join_vert(arma::zeros(n), arma::ones(2*m));

    const Vector xr = (result.solution.x.size() == n) ? result.solution.x : mu * arma::ones(n);

    result.solution.x = xr;

    ObjectiveFunction objective = [=](const Vector& x) mutable -> ObjectiveResult
    {
        auto xx = x.rows(0, n-1);
        auto xp = x.rows(n, n+m-1);
        auto xn = x.rows(n+m, n+2*m-1);
        objres.func = arma::sum(xp + xn) + 0.5*rho*arma::dot(xx - xr, xx - xr);
        objres.grad.rows(0, n-1) = rho*(xx - xr);
        return objres;
    };

    ConstraintFunction constraint = [=](const Vector& x) mutable -> ConstraintResult
    {
        auto xx = x.rows(0, n-1);
        auto xp = x.rows(n, n+m-1);
        auto xn = x.rows(n+m, n+2*m-1);

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
    result.solution.zl = options.ipnewton.mu/result.solution.x;

    ipnewton(problem, result, options);

    result.solution.x  = result.solution.x.rows(0, n-1);
    result.solution.y  = arma::zeros(m);
    result.solution.zl = result.solution.x/options.ipnewton.mu;
}

} // namespace Reaktor
