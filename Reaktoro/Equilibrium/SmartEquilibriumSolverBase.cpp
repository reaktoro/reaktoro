// Reaktoro is a unified framework for modeling chemically reactive systems.
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

#include <Reaktoro/Equilibrium/SmartEquilibriumSolverBase.hpp>

namespace Reaktoro {

/// Construct an SmartEquilibriumSolverBase instance with given chemical system.
SmartEquilibriumSolverBase::SmartEquilibriumSolverBase(const ChemicalSystem& system)
        : system(system), partition(Partition(system)),
          properties(system), solver(system)
{
    // Initialize the canonicalizer with the formula matrix Ae of the equilibrium species
    canonicalizer.compute(partition.formulaMatrixEquilibriumPartition());
}

/// Construct an SmartEquilibriumSolverBase instance with given partition of the chemical system.
SmartEquilibriumSolverBase::SmartEquilibriumSolverBase(const Partition& partition)
        : system(partition.system()), partition(partition),
          properties(partition.system()), solver(partition)
{
    // Initialize the canonicalizer with the formula matrix Ae of the equilibrium species
    canonicalizer.compute(partition.formulaMatrixEquilibriumPartition());

    // Initialize indices of the equilibrium and kinetic species
    ies = partition.indicesEquilibriumSpecies();
    iee = partition.indicesEquilibriumElements();
}

SmartEquilibriumSolverBase::~SmartEquilibriumSolverBase()
{
}

/// Set the options for the equilibrium calculation.
auto SmartEquilibriumSolverBase::setOptions(const SmartEquilibriumOptions& options_) -> void
{
    options = options_;

    // Tweak the options for the Gibbs energy minimization during learning operations.
    options.learning.hessian = GibbsHessian::Exact; // ensure the use of an exact Hessian of the Gibbs energy function
    options.learning.optimum.tolerance = 1e-10; // ensure the use of a stricter residual tolerance for the Gibbs energy minimization

    solver.setOptions(options.learning);
}

/// Solve the equilibrium problem with given problem definition
auto SmartEquilibriumSolverBase::solve(ChemicalState& state, const EquilibriumProblem& problem) -> SmartEquilibriumResult
{
    const auto T = problem.temperature();
    const auto P = problem.pressure();
    const auto& iee = partition.indicesEquilibriumElements();
    be = problem.elementAmounts()(iee);
    return solve(state, T, P, be);
}

auto SmartEquilibriumSolverBase::solve(ChemicalState& state, double T, double P, VectorConstRef be) -> SmartEquilibriumResult
{
    tic(SOLVE_STEP);

    // Absolutely ensure an exact Hessian of the Gibbs energy function is used in the calculations
    setOptions(options);

    // Reset the result of the last smart equilibrium calculation
    result = {};

    // Perform a smart estimate of the chemical state
    timeit( estimate(state, T, P, be), result.timing.estimate= );

    // Perform a learning step if the smart prediction is not sactisfatory
    if(!result.estimate.accepted)
        timeit( learn(state, T, P, be), result.timing.learn= );

    result.timing.solve = toc(SOLVE_STEP);

    return result;
}

auto SmartEquilibriumSolverBase::getProperties() const -> const ChemicalProperties&
{
    return properties;
}

auto SmartEquilibriumSolverBase::getResult() const -> const SmartEquilibriumResult&
{
    return result;
}

} // namespace Reaktoro

