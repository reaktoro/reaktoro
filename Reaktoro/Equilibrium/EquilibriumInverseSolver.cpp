// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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

#include "EquilibriumInverseSolver.hpp"

// Reaktoro includes
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Equilibrium/EquilibriumInverseProblem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktoro/Math/Roots.hpp>

namespace Reaktoro {

struct EquilibriumInverseSolver::Impl
{
    /// The chemical system instance
    ChemicalSystem system;

    /// The partition of the chemical system
    Partition partition;

    /// The options of the equilibrium solver
    EquilibriumOptions options;

    /// The solver for the equilibrium calculations
    EquilibriumSolver solver;

    /// Construct a Impl instance
    Impl(const ChemicalSystem& system)
    : system(system), solver(system)
    {
        // Set the default partition as all species are in equilibrium
        setPartition(Partition(system));
    }

    /// Set the partition of the chemical system
    auto setPartition(const Partition& partition_) -> void
    {
        // Set the partition of the chemical system
        partition = partition_;

        // Set the partition of the equilibrium solver
        solver.setPartition(partition);
    }

    /// Set the partition of the chemical system using a formatted string
    auto setPartition(std::string partition) -> void
    {
        setPartition(Partition(system, partition));
    }

    /// Solve an inverse equilibrium problem
    auto solve(ChemicalState& state, const EquilibriumInverseProblem& problem) -> EquilibriumResult
    {
        // The accumulated equilibrium result of this inverse problem calculation
        EquilibriumResult result;

        // Define auxiliary variables from the inverse problem definition
        const Matrix C = problem.formulaMatrixTitrants();
        const Vector b0 = problem.initialElementAmounts();
        const Indices ies = partition.indicesEquilibriumSpecies();
        const Indices iee = partition.indicesEquilibriumElements();

        // Get the rows corresponding to equilibrium elements only
        const Matrix Ce = rows(C, iee);
        const Vector be0 = rows(b0, iee);

        const Index num_titrants = C.cols();

        // Define auxiliary options
        const double tol = options.optimum.tolerance;
        const double maxiter = options.optimum.max_iterations;

        // Define auxiliary instances to avoid memory reallocation
        ChemicalProperties properties;
        ChemicalVector res;
        Matrix dndb;
        Matrix dfdn;

        // Set the options and partition in the equilibrium solver
        solver.setOptions(options);
        solver.setPartition(partition);

        // Define the non-linear residual function
        auto func = [&](const Vector& x, Vector& f, Matrix& J)
        {
            // The amounts of elements in the equilibrium partition
            const Vector be = be0 + Ce*x;

            // Solve the equilibrium problem with update `be`
            result += solver.solve(state, be);

            // Calculate the chemical properties of the system
            properties = state.properties();

            // Calculate the residuals of the equilibrium constraints
            res = problem.residualConstraints(properties);

            // Get the partial molar derivatives of f w.r.t. amounts of equilibrium species
            dfdn = cols(res.ddn, ies);

            // Set the residual and its derivatives w.r.t. x
            f = res.val;
            J = dfdn * dndb * C;
        };

        // Initialize the initial guess of the titrant amounts
        Vector x = zeros(num_titrants);

        // Apply Newton's method to find the amounts of titrants and
        // at the same time the equilibrium state of the system
        // satisfying the given activity and amount constraints
        x = newton(func, x, tol, maxiter);

        return result;
    }
};

EquilibriumInverseSolver::EquilibriumInverseSolver(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

EquilibriumInverseSolver::EquilibriumInverseSolver(const EquilibriumInverseSolver& other)
: pimpl(new Impl(*other.pimpl))
{}

EquilibriumInverseSolver::~EquilibriumInverseSolver()
{}

auto EquilibriumInverseSolver::operator=(EquilibriumInverseSolver other) -> EquilibriumInverseSolver&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto EquilibriumInverseSolver::setOptions(const EquilibriumOptions& options) -> void
{
    pimpl->options = options;
}

auto EquilibriumInverseSolver::setPartition(const Partition& partition) -> void
{
    pimpl->setPartition(partition);
}

auto EquilibriumInverseSolver::setPartition(std::string partition) -> void
{
    pimpl->setPartition(partition);
}

auto EquilibriumInverseSolver::solve(ChemicalState& state, const EquilibriumInverseProblem& problem) -> EquilibriumResult
{
    return pimpl->solve(state, problem);
}

} // namespace Reaktoro
