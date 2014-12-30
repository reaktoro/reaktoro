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

#include "EquilibriumSolver.hpp"

// Reaktor includes
#include <Reaktor/Core/ChemicalSystem.hpp>
#include <Reaktor/Core/Partition.hpp>
#include <Reaktor/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktor/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktor/Equilibrium/EquilibriumResult.hpp>
#include <Reaktor/Optimization/OptimumOptions.hpp>
#include <Reaktor/Optimization/OptimumProblem.hpp>
#include <Reaktor/Optimization/OptimumResult.hpp>
#include <Reaktor/Optimization/OptimumSolver.hpp>

namespace Reaktor {

struct EquilibriumSolver::Impl
{
    OptimumSolver optimum_solver;
};

EquilibriumSolver::EquilibriumSolver()
: pimpl(new Impl())
{}

EquilibriumSolver::EquilibriumSolver(const EquilibriumSolver& other)
: pimpl(new Impl(*other.pimpl))
{}

EquilibriumSolver::~EquilibriumSolver()
{}

auto EquilibriumSolver::operator=(EquilibriumSolver other) -> EquilibriumSolver&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto EquilibriumSolver::approximate(const EquilibriumProblem& problem, EquilibriumResult& result) -> void
{
    approximate(problem, result, {});
}

auto EquilibriumSolver::approximate(const EquilibriumProblem& problem, EquilibriumResult& result, const EquilibriumOptions& options) -> void
{
    OptimumProblem optimum_problem(problem);

    OptimumResult optimum_result;
    optimum_result.solution.x  = result.solution.n;
    optimum_result.solution.y  = result.solution.y;
    optimum_result.solution.zl = result.solution.z;

    pimpl->optimum_solver.approximate(optimum_problem, optimum_result, options.optimisation);

    result.solution.n = optimum_result.solution.x;
    result.solution.y = optimum_result.solution.y;
    result.solution.z = optimum_result.solution.zl;
    result.statistics = optimum_result.statistics;
}

auto EquilibriumSolver::solve(const EquilibriumProblem& problem, EquilibriumResult& result) -> void
{
    solve(problem, result, {});
}

auto EquilibriumSolver::solve(const EquilibriumProblem& problem, EquilibriumResult& result, const EquilibriumOptions& options) -> void
{
    OptimumProblem optimum_problem(problem);

    OptimumResult optimum_result;
    optimum_result.solution.x  = result.solution.n;
    optimum_result.solution.y  = result.solution.y;
    optimum_result.solution.zl = result.solution.z;

    OptimumOptions optimum_options = options.optimisation;

    if(options.hessian == DiagonalHessian)
    {
        optimum_options.ipnewton.kkt.algorithm = Reaktor::KktRangespace;
        optimum_options.ipopt.kkt.algorithm = Reaktor::KktRangespace;
        optimum_options.ipnewton.kkt.diagonalH = true;
        optimum_options.ipopt.kkt.diagonalH = true;
    }

    pimpl->optimum_solver.solve(optimum_problem, optimum_result, optimum_options);

    result.solution.n = optimum_result.solution.x;
    result.solution.y = optimum_result.solution.y;
    result.solution.z = optimum_result.solution.zl;
    result.statistics = optimum_result.statistics;
}

} // namespace Reaktor
