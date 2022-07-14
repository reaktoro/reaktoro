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

#include "SmartEquilibriumSolver.hpp"

// Reaktoro includes
#include <Reaktoro/Equilibrium/SmartEquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumSolverClustering.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumSolverPriorityQueue.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumSolverNN.hpp>

namespace Reaktoro {

SmartEquilibriumSolver::SmartEquilibriumSolver(const ChemicalSystem& system) {
    solverptr = std::make_unique<SmartEquilibriumSolverClustering>(system);
}

SmartEquilibriumSolver::SmartEquilibriumSolver(const Partition& partition) {
    solverptr = std::make_unique<SmartEquilibriumSolverClustering>(partition);
}

SmartEquilibriumSolver::~SmartEquilibriumSolver()
{
}

/// Set the options for the equilibrium calculation.
auto SmartEquilibriumSolver::setOptions(const SmartEquilibriumOptions& options) -> void
{
    if(!initialized)
    {
        // Switch the tag as smart equilibrium solver is initialized
        initialized = true;

        // Initialize smart equilibrium solver depending on the chosen options
        switch(options.method)
        {
            case SmartEquilibriumStrategy::NearestNeighbour:
                solverptr = std::make_unique<SmartEquilibriumSolverNN>(solverptr->partition());
                break;
            case SmartEquilibriumStrategy::PriorityQueue:
                solverptr = std::make_unique<SmartEquilibriumSolverPriorityQueue>(solverptr->partition());
                break;
            case SmartEquilibriumStrategy::Clustering:
                solverptr = std::make_unique<SmartEquilibriumSolverClustering>(solverptr->partition());
                break;
            default:
                solverptr = std::make_unique<SmartEquilibriumSolverClustering>(solverptr->partition());
                break;
        }
    }

    // Set options of the smart equilibrium options
    solverptr->setOptions(options);
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, const EquilibriumProblem& problem) -> SmartEquilibriumResult
{
    return solverptr->solve(state, problem);
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, double T, double P, VectorConstRef be) -> SmartEquilibriumResult
{
    return solverptr->solve(state, T, P, be);
}

auto SmartEquilibriumSolver::result() const -> const SmartEquilibriumResult&
{
    return solverptr->result();
}

auto SmartEquilibriumSolver::outputInfo() const -> void
{
    solverptr->outputInfo();
}

} // namespace Reaktoro
