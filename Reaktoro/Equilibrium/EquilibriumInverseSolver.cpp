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

    /// The options of the equilibrium solver
    EquilibriumOptions options;

    /// The solver for the equilibrium calculations
    EquilibriumSolver solver;

    /// Construct a Impl instance
    Impl(const ChemicalSystem& system)
    : system(system), solver(system)
    {
    }

    /// Solve an inverse equilibrium problem
    auto solve(ChemicalState& state, const EquilibriumInverseProblem& problem) -> EquilibriumResult
    {
        // The accumulated equilibrium result of this inverse problem calculation
        EquilibriumResult result;

        // Define auxiliary variables from the inverse problem definition
        const Matrix C = problem.formulaMatrixTitrants();
        const Vector b0 = problem.initialElementAmounts();
        const Indices ainds = problem.indicesSpeciesActivityConstraints();
        const Indices ninds = problem.indicesSpeciesAmountConstraints();
        const Vector avals = problem.valuesSpeciesActivityConstraints();
        const Vector nvals = problem.valuesSpeciesAmountConstraints();
        const Partition partition = problem.partition();
        const Index num_titrants = C.cols();

        // Define auxiliary options
        const double tol = options.optimum.tolerance;
        const double maxiter = options.optimum.max_iterations;

        // Define auxiliary instances to avoid memory reallocation
        ChemicalProperties properties;
        ChemicalVector a;
        Matrix dndb;
        Vector n;

        // Set the options and partition in the equilibrium solver
        solver.setOptions(options);
        solver.setPartition(partition);

        // Define the non-linear residual function
        auto func = [&](const Vector& x, Vector& f, Matrix& J)
        {
            const Vector b = b0 + C*x;

            result += solver.solve(state, b);

            properties = state.properties();

            n = state.speciesAmounts();
            a = exp(properties.lnActivities());

            Index offset = 0;

            // Loop over all species activity constraints
            for(Index i = 0; i < ainds.size(); ++i)
            {
                // Get the index of current species with fixed activity
                const Index ispecies = ainds[i];

                // Compute the residual of the activity constraint
                f[offset + i] = a[ispecies].val - avals[i];

                // Compute the x-derivatives of the residual of the activity constraint
                J.row(offset + i) = tr(a[ispecies].ddn) * dndb * C;
            }

            // Update the offset variable
            offset += ainds.size();

            // Loop over all species amount constraints
            for(Index i = 0; i < ninds.size(); ++i)
            {
                // Get the index of current species with fixed amount
                const Index ispecies = ninds[i];

                // Compute the residual of the amount constraint
                f[offset + i] = n[ispecies] - nvals[i];

                // Compute the x-derivatives of the residual of the amount constraint
                J.row(offset + i) = tr(dndb.row(ispecies)) * C;
            }
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

auto EquilibriumInverseSolver::solve(ChemicalState& state, const EquilibriumInverseProblem& problem) -> EquilibriumResult
{
    return pimpl->solve(state, problem);
}

} // namespace Reaktoro
