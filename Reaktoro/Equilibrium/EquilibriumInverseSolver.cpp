// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#include "EquilibriumInverseSolver.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Equilibrium/EquilibriumInverseProblem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktoro/Optimization/NonlinearSolver.hpp>

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
    EquilibriumSolver equilibriumsolver;

    /// The sensitivity of the equilibrium state
    EquilibriumSensitivity sensitivity;

    /// Construct a Impl instance with given chemical system
    Impl(const ChemicalSystem& system)
    : Impl(Partition(system))
    {
    }

    /// Construct a Impl instance with given partition of chemical system
    Impl(const Partition& partition)
    : system(partition.system()),
      partition(partition),
      equilibriumsolver(partition)
    {
    }

    /// Solve an inverse equilibrium problem
    auto solve(ChemicalState& state, const EquilibriumInverseProblem& problem) -> EquilibriumResult
    {
        // The accumulated equilibrium result of this inverse problem calculation
        EquilibriumResult result;

        // Define auxiliary variables from the inverse problem definition
        const Index Nt = problem.numTitrants();
        const Index Nc = problem.numConstraints();
        const Matrix C = problem.formulaMatrixTitrants();
        const Vector b0 = problem.elementInitialAmounts();
        const Indices ies = partition.indicesEquilibriumSpecies();
        const Indices iee = partition.indicesEquilibriumElements();

        // Get the rows corresponding to equilibrium elements only
        const Matrix Ce = rows(C, iee);
        const Vector be0 = rows(b0, iee);

        // The temperature and pressure for the calculation
        const double T = problem.temperature();
        const double P = problem.pressure();

        // Set the temperature and pressure of the chemical state
        state.setTemperature(T);
        state.setPressure(P);

        // Define auxiliary instances to avoid memory reallocation
        ChemicalProperties properties;
        ResidualEquilibriumConstraints res;
        NonlinearResidual nonlinear_residual;

        // Auxiliary references to the non-linear residual data
        auto& F = nonlinear_residual.val;
        auto& J = nonlinear_residual.jacobian;

        // Set the options and partition in the equilibrium solver
        equilibriumsolver.setOptions(options);

        // Define the non-linear problem with inequality constraints
        NonlinearProblem nonlinear_problem;

        // Set the linear inequality constraints of the titrant molar amounts
        nonlinear_problem.n = Nt;
        nonlinear_problem.m = Nc;
        nonlinear_problem.A = C;
        nonlinear_problem.b = -be0;

        // Set the non-linear function of the non-linear problem
        nonlinear_problem.f = [&](VectorConstRef x) mutable
        {
            // The amounts of elements in the equilibrium partition
            const Vector be = be0 + Ce*x;

            // Solve the equilibrium problem with update `be`
            result += equilibriumsolver.solve(state, T, P, be);

            // Check if the equilibrium calculation converged
            if(!result.optimum.succeeded)
            {
                // If not, solve using cold start
                state.setSpeciesAmounts(0.0);
                result += equilibriumsolver.solve(state, T, P, be);
            }

            // Check if the function evaluation was successful
            nonlinear_residual.succeeded = result.optimum.succeeded;

            // Update the sensitivity of the equilibrium state
            sensitivity = equilibriumsolver.sensitivity();

            // Calculate the residuals of the equilibrium constraints
            res = problem.residualEquilibriumConstraints(x, state);

            // Calculate the residual vector `F` and its Jacobian `J`
            F = res.val;
            J = res.ddx + res.ddn * sensitivity.dndb * C;

            return nonlinear_residual;
        };

        // Initialize the initial guess of the titrant amounts
        Vector x = problem.titrantInitialAmounts();

        // Replace zeros in x by small molar amounts
        x = (x.array() > 0.0).select(x, 1e-6);

        // Solve the non-linear problem with inequality constraints
        NonlinearSolver nonlinear_solver;
        nonlinear_solver.solve(nonlinear_problem, x, options.nonlinear);

        return result;
    }
};

EquilibriumInverseSolver::EquilibriumInverseSolver(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

EquilibriumInverseSolver::EquilibriumInverseSolver(const Partition& partition)
: pimpl(new Impl(partition))
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
    RuntimeError("Cannot proceed with EquilibriumInverseSolver::setPartition.",
        "EquilibriumInverseSolver::setPartition is deprecated. "
        "Use constructor EquilibriumInverseSolver(const Partition&) instead.");
}

auto EquilibriumInverseSolver::solve(ChemicalState& state, const EquilibriumInverseProblem& problem) -> EquilibriumResult
{
    return pimpl->solve(state, problem);
}

auto EquilibriumInverseSolver::sensitivity() -> EquilibriumSensitivity
{
    return pimpl->equilibriumsolver.sensitivity();
}

} // namespace Reaktoro
