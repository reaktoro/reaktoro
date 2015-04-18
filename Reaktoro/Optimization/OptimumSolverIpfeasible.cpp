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

#include "OptimumSolverIpfeasible.hpp"

// Reaktoro includes
#include <Reaktoro/Optimization/OptimumOptions.hpp>
#include <Reaktoro/Optimization/OptimumProblem.hpp>
#include <Reaktoro/Optimization/OptimumResult.hpp>
#include <Reaktoro/Optimization/OptimumState.hpp>
#include <Reaktoro/Optimization/OptimumSolverIpnewton.hpp>

namespace Reaktoro {

struct OptimumSolverIpfeasible::Impl
{
    OptimumSolverIpnewton ipnewton;

    auto approximate(OptimumProblem problem, OptimumState& state, OptimumOptions options) -> OptimumResult;
};

auto OptimumSolverIpfeasible::Impl::approximate(OptimumProblem problem, OptimumState& state, OptimumOptions options) -> OptimumResult
{
    // Auxiliary variables
    const unsigned n = problem.A.cols();
    const unsigned m = problem.A.rows();
    const double  mu = options.ipnewton.mu;
    const double rho = std::sqrt(mu);

    // The total number of primal variables for the feasibility problem
    const unsigned t = n + 2*m;

    // Initialize the solution variables if it has not been done before
    if(state.x.size() != n)
        state.x = mu * ones(n);

    // The reference point from which the solution cannot differ much
    const Vector xr = state.x;

    // The definition of the feasibility problem
    OptimumProblem fproblem;

    // The result of the objective function evaluation of the feasibility problem
    ObjectiveResult res;

    // Initialize the gradient member
    res.grad.resize(t);
    res.grad << zeros(n), ones(2*m);

    // Initialize the hessian member
    res.hessian.mode = Hessian::Diagonal;
    res.hessian.diagonal = zeros(t);
    rows(res.hessian.diagonal, 0, n) = rho * ones(n);

    // Define the objective function of the feasibility problem
    fproblem.objective = [=](const Vector& x) mutable
    {
        const auto xx = rows(x, 0, n);
        const auto xp = rows(x, n, m);
        const auto xn = rows(x, n + m, m);
        res.val = (xp + xn).sum() + 0.5 * rho * (xx - xr).dot(xx - xr);
        rows(res.grad, 0, n) = rho*(xx - xr);
        return res;
    };

    // Define the equality constraint of the feasibility problem
    fproblem.l = zeros(t);
    fproblem.b = problem.b;
    fproblem.A.resize(m, t);
    cols(fproblem.A, 0, n)     =  problem.A;
    cols(fproblem.A, n, m)     = -identity(m, m);
    cols(fproblem.A, n + m, m) =  identity(m, m);

    // Set the initial guess
    const Vector xx = state.x;
    const Vector xp = mu * ones(m);
    const Vector xn = mu * ones(m);

    state.x.resize(t);
    state.x << xx, xp, xn;
    state.y  = zeros(m);
    state.z = mu/state.x.array();

    // Solve the feasibility problem
    auto result = ipnewton.solve(fproblem, state, options);

    // Convert the result of the artificial feasibility problem
    state.x = rows(state.x, 0, n);
    state.y = zeros(m);
    state.z = mu/state.x.array();

    return result;
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

auto OptimumSolverIpfeasible::approximate(const OptimumProblem& problem, OptimumState& state) -> OptimumResult
{
    return approximate(problem, state, {});
}

auto OptimumSolverIpfeasible::approximate(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
{
    return pimpl->approximate(problem, state, options);
}

} // namespace Reaktoro
