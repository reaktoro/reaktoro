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

#pragma once

// C++ includes
#include <memory>

// Reaktoro includes
#include <Reaktoro/Equilibrium/SmartEquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumSolverBase.hpp>

namespace Reaktoro {

/// A class used to perform equilibrium calculations using machine learning scheme.
class SmartEquilibriumSolver
{
public:
    /// Construct an SmartEquilibriumSolverBase instance with given chemical system
    explicit SmartEquilibriumSolver(const ChemicalSystem& system);

    /// Construct an SmartEquilibriumSolverBase instance with given partition
    explicit SmartEquilibriumSolver(const Partition& partition);

    /// Destroy this SmartEquilibriumSolverBase instance.
    virtual ~SmartEquilibriumSolver();

    /// Set the options for the equilibrium calculation.
    auto setOptions(const SmartEquilibriumOptions& options) -> void;

    /// Solve a chemical equilibrium problem.
    /// @param state[in,out] The initial guess and the final state of the equilibrium calculation
    /// @param T The temperature (in units of K)
    /// @param P The pressure (in units of Pa)
    /// @param be The amounts of the elements in the equilibrium partition
    auto solve(ChemicalState& state, double T, double P, VectorConstRef be) -> SmartEquilibriumResult;

    /// Solve a chemical equilibrium problem.
    /// @param state[in,out] The initial guess and the final state of the equilibrium calculation
    /// @param problem The equilibrium problem with given temperature, pressure, and element amounts.
    auto solve(ChemicalState& state, const EquilibriumProblem& problem) -> SmartEquilibriumResult;

    /// Return the result of the last smart equilibrium calculation.
    auto result() const -> const SmartEquilibriumResult&;

    /// Output into about the ODML algorithm
    auto outputInfo() const -> void;

private:
    // Pointer to the smart equilibrium solver instance
    std::unique_ptr<SmartEquilibriumSolverBase> solverptr;

    // Flag to indicate whether the solver has been initialized
    bool initialized = false;
};

} // namespace Reaktoro
