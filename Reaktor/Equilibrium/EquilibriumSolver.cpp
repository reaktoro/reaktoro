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

auto convert(const EquilibriumProblem& problem, EquilibriumResult& result) -> OptimumProblem
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
    auto update = [=](EquilibriumResult& result, const Vector& x)
    {
        // Auxiliary references to the members of `result`
        Vector& n  = result.solution.n;
        Vector& u  = result.solution.u;
        Vector& ne = result.solution.ne;
        Vector& ue = result.solution.ue;
        double& ge = result.solution.ge;

        // Set the molar amounts of the species
        rows(n, iequilibrium) = x;

        // Set the molar amounts of the equilibrium species
        ne.noalias() = x;

        // Set the scaled chemical potentials of the species
        u = system.chemicalPotentials(T, P, n).val()/RT;

        // Set the scaled chemical potentials of the equilibrium species
        rows(u, iequilibrium).to(ue);

        // Set the Gibbs energy of the equilibrium partition
        ge = dot(ne, ue);
    };

    // Define the Gibbs energy function
    ObjectiveFunction gibbs = [=,&result](const Vector& x) mutable
    {
        update(result, x);
        return result.solution.ge;
    };

    // Define the gradient of the Gibbs energy function
    ObjectiveGradFunction gibbs_grad = [=,&result](const Vector& x) mutable
    {
        if(x == result.solution.ne)
            return result.solution.ue;
        update(result, x);
        return result.solution.ue;
    };

    // Define the Hessian of the Gibbs energy function
    ObjectiveDiagonalHessianFunction gibbs_hessian = [=](const Vector& n) mutable
    {
        return inv(n);
    };

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
    OptimumProblem optimum_problem = convert(problem, result);

    OptimumResult optimum_result;
    optimum_result.solution.x  = result.solution.n;
    optimum_result.solution.y  = result.solution.ye;
    optimum_result.solution.z = result.solution.ze;

    pimpl->optimum_solver.approximate(optimum_problem, optimum_result, options.optimum);

    result.solution.ne = optimum_result.solution.x;
    result.solution.ye = optimum_result.solution.y;
    result.solution.ze = optimum_result.solution.z;
    result.statistics = optimum_result.statistics;
}

auto EquilibriumSolver::solve(const EquilibriumProblem& problem, EquilibriumResult& result) -> void
{
    solve(problem, result, {});
}

auto EquilibriumSolver::solve(const EquilibriumProblem& problem, EquilibriumResult& result, const EquilibriumOptions& options) -> void
{
    // The reference to the chemical system and its partitioning
    const Partition& partition = problem.partition();

    // The indices of the equilibrium species
    const Indices& iequilibrium = partition.equilibriumSpeciesIndices();

    // Ensure the initial guess of the primal variables (i.e., the molar amounts
    // of the equilibrium species) is extracted from the molar amounts of the species
    result.optimum.solution.x = rows(result.solution.n, iequilibrium);

    // Convert an EquilibriumProblem into an OptimumProblem
    OptimumProblem optimum_problem = convert(problem, result);

    // Solve the OptimumProblem
    pimpl->optimum_solver.solve(optimum_problem, result.optimum, options.optimum);

    // Copy the statistics of the optimisation calculation
    result.statistics = result.optimum.statistics;
}

} // namespace Reaktor
