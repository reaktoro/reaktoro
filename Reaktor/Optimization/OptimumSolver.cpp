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

#include "OptimumSolver.hpp"

// Reaktor includes
#include <Reaktor/Optimization/OptimumOptions.hpp>
#include <Reaktor/Optimization/OptimumProblem.hpp>
#include <Reaktor/Optimization/OptimumResult.hpp>
#include <Reaktor/Optimization/OptimumSolverIpfeasible.hpp>
#include <Reaktor/Optimization/OptimumSolverIpnewton.hpp>
#include <Reaktor/Optimization/OptimumSolverIpopt.hpp>

namespace Reaktor {

struct OptimumSolver::Impl
{
    OptimumSolverIpfeasible ipfeasible;

    OptimumSolverIpnewton ipnewton;

    OptimumSolverIpopt ipopt;
};

OptimumSolver::OptimumSolver()
: pimpl(new Impl())
{}

OptimumSolver::OptimumSolver(const OptimumSolver& other)
: pimpl(new Impl(*other.pimpl))
{}

OptimumSolver::~OptimumSolver()
{}

auto OptimumSolver::operator=(OptimumSolver other) -> OptimumSolver&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto OptimumSolver::approximate(const OptimumProblem& problem, OptimumState& state) -> OptimumResult
{
    return approximate(problem, state, {});
}

auto OptimumSolver::approximate(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
{
    return pimpl->ipfeasible.approximate(problem, state, options);
}

auto OptimumSolver::solve(const OptimumProblem& problem, OptimumState& state) -> OptimumResult
{
    return solve(problem, state, {});
}

auto OptimumSolver::solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
{
    switch(options.method)
    {
    case OptimumMethod::Ipnewton:
        return pimpl->ipnewton.solve(problem, state, options);
    case OptimumMethod::Ipopt:
        return pimpl->ipopt.solve(problem, state, options);
    default:
        return pimpl->ipnewton.solve(problem, state, options);
    }
}

} // namespace Reaktor
