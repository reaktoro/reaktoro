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
#include <Reaktor/Common/Constants.hpp>
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
namespace {

auto convert(const EquilibriumProblem& problem) -> OptimumProblem
{
    // The reference to the chemical system and its partitioning
    const ChemicalSystem& system = problem.system();
    const Partition& partition = problem.partition();

    // The number of equilibrium species and linearly independent components in the equilibrium partition
    const unsigned num_equilibrium_species = partition.equilibriumSpeciesIndices().size();
    const unsigned num_components = problem.components().size();

    // The temperature and pressure of the equilibrium calculation
    const double T  = problem.temperature();
    const double P  = problem.pressure();
    const double RT = universalGasConstant*T;

    // An auxiliary vector to hold the chemical potentials of the species
    // scaled by RT. This is a shared pointer because (1) it must outlive
    // this function and (2) its content is shared among the functions below
    std::shared_ptr<Vector> u_ptr(new Vector());

    // Define the Gibbs energy function
    ObjectiveFunction gibbs = [=](const Vector& n) mutable
    {
        *u_ptr = (1/RT)*system.chemicalPotentials(T, P, n).val();
        return dot(n, *u_ptr);
    };

    // Define the gradient of the Gibbs energy function
    ObjectiveGradFunction gibbs_grad = [=](const Vector& n) mutable
    {
        return *u_ptr;
    };

    // Define the Hessian of the Gibbs energy function
    ObjectiveDiagonalHessianFunction gibbs_hessian = [=](const Vector& n) mutable
    {
        return inv(n);
    };

    // The left-hand side matrix of the linearly independent mass-charge balance equations
    const Matrix A = problem.balanceMatrix();

    // The right-hand side vector of the linearly independent mass-charge balance equations
    const Vector b = problem.componentAmounts();

    // Define the mass-cahrge balance contraint function
    ConstraintFunction balance_constraint = [=](const Vector& n) -> Vector
    {
        return A*n - b;
    };

    // Define the gradient function of the mass-cahrge balance contraint function
    ConstraintGradFunction balance_constraint_grad = [=](const Vector& n) -> Matrix
    {
        return A;
    };

    // Setup an OptimumProblem instance with the Gibbs energy function and the balance constraints
    OptimumProblem optimum_problem(num_equilibrium_species, num_components);
    optimum_problem.setObjective(gibbs);
    optimum_problem.setObjectiveGrad(gibbs_grad);
    optimum_problem.setObjectiveDiagonalHessian(gibbs_hessian);
    optimum_problem.setConstraint(balance_constraint);
    optimum_problem.setConstraintGrad(balance_constraint_grad);
    optimum_problem.setLowerBounds(0.0);

    return optimum_problem;
}

} // namespace

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
    OptimumProblem optimum_problem = convert(problem);

    OptimumResult optimum_result;
    optimum_result.solution.x  = result.solution.n;
    optimum_result.solution.y  = result.solution.y;
    optimum_result.solution.z = result.solution.z;

    pimpl->optimum_solver.approximate(optimum_problem, optimum_result, options.optimisation);

    result.solution.n = optimum_result.solution.x;
    result.solution.y = optimum_result.solution.y;
    result.solution.z = optimum_result.solution.z;
    result.statistics = optimum_result.statistics;
}

auto EquilibriumSolver::solve(const EquilibriumProblem& problem, EquilibriumResult& result) -> void
{
    solve(problem, result, {});
}

auto EquilibriumSolver::solve(const EquilibriumProblem& problem, EquilibriumResult& result, const EquilibriumOptions& options) -> void
{
    OptimumProblem optimum_problem = convert(problem);

    OptimumResult optimum_result;
    optimum_result.solution.x  = result.solution.n;
    optimum_result.solution.y  = result.solution.y;
    optimum_result.solution.z = result.solution.z;

    OptimumOptions optimum_options = options.optimisation;

    pimpl->optimum_solver.solve(optimum_problem, optimum_result, optimum_options);

    result.solution.n = optimum_result.solution.x;
    result.solution.y = optimum_result.solution.y;
    result.solution.z = optimum_result.solution.z;
    result.statistics = optimum_result.statistics;
}

} // namespace Reaktor
