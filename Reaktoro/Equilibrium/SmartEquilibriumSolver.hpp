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
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

// Forward declarations (classes)
class ChemicalProperties;
class ChemicalState;
class ChemicalSystem;
class EquilibriumProblem;
class EquilibriumSensitivity;
class Partition;

// Forward declarations (structs)
struct SmartEquilibriumOptions;
struct SmartEquilibriumResult;

/// A class used to perform equilibrium calculations using machine learning scheme.
class SmartEquilibriumSolver
{
public:
    /// Construct a default SmartEquilibriumSolver instance.
    SmartEquilibriumSolver();

    /// Construct an SmartEquilibriumSolver instance with given chemical system
    explicit SmartEquilibriumSolver(const ChemicalSystem& system);

    /// Construct an SmartEquilibriumSolver instance with given partition
    explicit SmartEquilibriumSolver(const Partition& partition);

    /// Construct a copy of an SmartEquilibriumSolver instance.
    SmartEquilibriumSolver(const SmartEquilibriumSolver& other);

    /// Assign an SmartEquilibriumSolver instance to this.
    auto operator=(SmartEquilibriumSolver other) -> SmartEquilibriumSolver&;

    /// Destroy this SmartEquilibriumSolver instance.
    virtual ~SmartEquilibriumSolver();

    /// Set the options for the equilibrium calculation.
    auto setOptions(const SmartEquilibriumOptions& options) -> void;

    /// Set the partition of the chemical system.
    auto setPartition(const Partition& partition) -> void;

    /// Solve a chemical equilibrium problem.
    /// @param state[in,out] The initial guess and the final state of the equilibrium calculation
    /// @param T The temperature (in units of K)
    /// @param P The pressure (in units of Pa)
    /// @param be The amounts of the elements in the equilibrium partition
    auto solve(ChemicalState& state, double T, double P, VectorConstRef be) -> SmartEquilibriumResult;

    auto solve(ChemicalState& state, double T, double P, VectorConstRef be, Index istep, Index icell) -> SmartEquilibriumResult;

    /// Solve a chemical equilibrium problem.
    /// @param state[in,out] The initial guess and the final state of the equilibrium calculation
    /// @param problem The equilibrium problem with given temperature, pressure, and element amounts.
    auto solve(ChemicalState& state, const EquilibriumProblem& problem) -> SmartEquilibriumResult;

    /// Return the chemical properties of the calculated equilibrium state.
    auto properties() const -> const ChemicalProperties&;

    /// Return the result of the last smart equilibrium calculation.
    auto result() const -> const SmartEquilibriumResult&;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
