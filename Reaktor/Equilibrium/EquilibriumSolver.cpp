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

struct EquilibriumSolver::Impl
{
    /// The solver for the optimisation calculations
    OptimumSolver optimum_solver;

    /// The result of an optimisation calculation
    OptimumResult optimum_result;

    /// The molar amounts of the equilibrium species
    Vector ne;

    /// The chemical potentials of the species
    Vector u;

    /// The chemical potentials of the equilibrium species
    Vector ue;

    /// Convert an EquilibriumProblem into an OptimumProblem
    auto convert(const EquilibriumProblem& problem, EquilibriumResult& result) -> OptimumProblem;

    /// Find a initial guess for an equilibrium problem
    auto approximate(const EquilibriumProblem& problem, EquilibriumResult& result, const EquilibriumOptions& options) -> void;

    /// Solve the equilibrium problem
    auto solve(const EquilibriumProblem& problem, EquilibriumResult& result, const EquilibriumOptions& options) -> void;
};

auto EquilibriumSolver::Impl::convert(const EquilibriumProblem& problem, EquilibriumResult& result) -> OptimumProblem
{
    // The reference to the chemical system and its partitioning
    const ChemicalSystem& system = problem.system();
    const Partition& partition = problem.partition();

    // The indices of the equilibrium species
    const Indices& iequilibrium = partition.equilibriumSpeciesIndices();

    // The number of equilibrium species and linearly independent components in the equilibrium partition
    const unsigned num_equilibrium_species = iequilibrium.size();
    const unsigned num_components = problem.components().size();

    // The temperature and pressure of the equilibrium calculation
    const double T  = problem.temperature();
    const double P  = problem.pressure();
    const double RT = universalGasConstant*T;

    // The left-hand side matrix of the linearly independent mass-charge balance equations
    const Matrix A = problem.balanceMatrix();

    // The right-hand side vector of the linearly independent mass-charge balance equations
    const Vector b = problem.componentAmounts();

    // Auxiliary function to update the state of a EquilibriumResult instance
    auto update = [=](EquilibriumResult& result, const Vector& x) mutable
    {
        // Set the molar amounts of the species
        rows(result.n, iequilibrium) = x;

        // Set the molar amounts of the equilibrium species
        ne.noalias() = x;

        // Set the scaled chemical potentials of the species
        u = system.chemicalPotentials(T, P, result.n).val()/RT;

        // Set the scaled chemical potentials of the equilibrium species
        rows(u, iequilibrium).to(ue);
    };

    // Define the Gibbs energy function
    ObjectiveFunction gibbs = [=,&result](const Vector& x) mutable
    {
        update(result, x);
        return dot(ne, ue);
    };

    // Define the gradient of the Gibbs energy function
    ObjectiveGradFunction gibbs_grad = [=,&result](const Vector& x) mutable
    {
        if(x == ne) return ue;
        update(result, x);
        return ue;
    };

    // Define the Hessian of the Gibbs energy function
    ObjectiveDiagonalHessianFunction gibbs_hessian = [=](const Vector& x) mutable
    {
        return inv(x);
    };

    // Define the mass-cahrge balance contraint function
    ConstraintFunction balance_constraint = [=](const Vector& x) -> Vector
    {
        return A*x - b;
    };

    // Define the gradient function of the mass-cahrge balance contraint function
    ConstraintGradFunction balance_constraint_grad = [=](const Vector& x) -> Matrix
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

auto EquilibriumSolver::Impl::approximate(const EquilibriumProblem& problem, EquilibriumResult& result, const EquilibriumOptions& options) -> void
{
    // Convert an EquilibriumProblem into an OptimumProblem
    OptimumProblem optimum_problem = convert(problem, result);

    // The reference to the chemical system and its partitioning
    const Partition& partition = problem.partition();

    // The indices of the equilibrium species
    const Indices& iequilibrium = partition.equilibriumSpeciesIndices();

    // Ensure the initial guess of the primal variables (i.e., the molar amounts
    // of the equilibrium species) is extracted from the molar amounts of the species
    rows(result.n, iequilibrium).to(result.optimum.solution.x);

    // Find an approximation to the optimisation problem
    optimum_solver.approximate(optimum_problem, optimum_result, options.optimum);

    // Copy the optimisation result to the equilibrium result
    rows(result.n, iequilibrium) = result.optimum.solution.x;

    // Copy the statistics of the optimisation calculation
    result.statistics = result.optimum.statistics;
}

auto EquilibriumSolver::Impl::solve(const EquilibriumProblem& problem, EquilibriumResult& result, const EquilibriumOptions& options) -> void
{
    // Convert an EquilibriumProblem into an OptimumProblem
    OptimumProblem optimum_problem = convert(problem, result);

    // The reference to the chemical system and its partitioning
    const Partition& partition = problem.partition();

    // The indices of the equilibrium species
    const Indices& iequilibrium = partition.equilibriumSpeciesIndices();

    // Ensure the initial guess of the primal variables (i.e., the molar amounts
    // of the equilibrium species) is extracted from the molar amounts of the species
    rows(result.n, iequilibrium).to(result.optimum.solution.x);

    // Solve the optimisation problem
    optimum_solver.solve(optimum_problem, result.optimum, options.optimum);

    // Copy the optimisation result to the equilibrium result
    rows(result.n, iequilibrium) = result.optimum.solution.x;

    // Copy the statistics of the optimisation calculation
    result.statistics = result.optimum.statistics;
}

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
    pimpl->approximate(problem, result, options);
}

auto EquilibriumSolver::solve(const EquilibriumProblem& problem, EquilibriumResult& result) -> void
{
    solve(problem, result, {});
}

auto EquilibriumSolver::solve(const EquilibriumProblem& problem, EquilibriumResult& result, const EquilibriumOptions& options) -> void
{
    pimpl->solve(problem, result, options);
}

} // namespace Reaktor
