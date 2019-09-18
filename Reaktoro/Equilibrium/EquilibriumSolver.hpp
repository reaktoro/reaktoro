// This file is part of Reaktoro (https://reaktoro.org).
//
// Reaktoro is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// Reaktoro is distributed in the hope that it will be useful,
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

// Forward declarations
class ChemicalProperties;
class ChemicalState;
class ChemicalSystem;
class Partition;
class EquilibriumProblem;
struct EquilibriumOptions;
struct EquilibriumResult;
struct EquilibriumSensitivity;

/// A solver class for solving chemical equilibrium calculations.
class EquilibriumSolver
{
public:
    /// Construct a default EquilibriumSolver instance
    EquilibriumSolver();

    /// Construct an EquilibriumSolver instance
    explicit EquilibriumSolver(const ChemicalSystem& system);

    /// Construct an EquilibriumSolver instance with given partition
    explicit EquilibriumSolver(const Partition& partition);

    /// Construct a copy of an EquilibriumSolver instance
    EquilibriumSolver(const EquilibriumSolver& other);

    /// Destroy this EquilibriumSolver instance
    virtual ~EquilibriumSolver();

    /// Assign a copy of an EquilibriumSolver instance
    auto operator=(EquilibriumSolver other) -> EquilibriumSolver&;

    /// Set the options of the equilibrium solver
    auto setOptions(const EquilibriumOptions& options) -> void;

    /// Set the partition of the chemical system
    auto setPartition(const Partition& partition) -> void;

    /// Find an initial feasible guess for an equilibrium problem.
    /// @param state[in,out] The initial guess and the final state of the equilibrium approximation
    /// @param be The molar amounts of the elements in the equilibrium partition
    auto approximate(ChemicalState& state, double T, double P, VectorConstRef be) -> EquilibriumResult;

    /// Find an initial feasible guess for an equilibrium problem.
    /// @param state[in,out] The initial guess and the final state of the equilibrium approximation
    /// @param problem The equilibrium problem with given temperature, pressure, and element amounts.
    auto approximate(ChemicalState& state, const EquilibriumProblem& problem) -> EquilibriumResult;

    /// Find an initial feasible guess for an equilibrium problem.
    /// @param state[in,out] The initial guess and the final state of the equilibrium approximation
    auto approximate(ChemicalState& state) -> EquilibriumResult;

    /// Solve an equilibrium problem with given molar amounts of the elements in the equilibrium partition..
    /// @param state[in,out] The initial guess and the final state of the equilibrium calculation
    /// @param be The molar amounts of the elements in the equilibrium partition
    auto solve(ChemicalState& state, double T, double P, VectorConstRef be) -> EquilibriumResult;

    /// Solve an equilibrium problem with given molar amounts of the elements in the equilibrium partition..
    /// @param state[in,out] The initial guess and the final state of the equilibrium calculation
    /// @param T The temperature (in units of K)
    /// @param P The pressure (in units of Pa)
    /// @param be The molar amounts of the elements in the equilibrium partition
    auto solve(ChemicalState& state, double T, double P, const double* be) -> EquilibriumResult;

    /// Solve an equilibrium problem with given equilibrium problem.
    /// @param state[in,out] The initial guess and the final state of the equilibrium calculation
    /// @param problem The equilibrium problem with given temperature, pressure, and element amounts.
    auto solve(ChemicalState& state, const EquilibriumProblem& problem) -> EquilibriumResult;

    /// Solve an equilibrium problem with given chemical state in disequilibrium.
    /// @param state[in,out] The initial guess and the final state of the equilibrium calculation
    auto solve(ChemicalState& state) -> EquilibriumResult;

    /// Return the chemical properties of the calculated equilibrium state.
    auto properties() const -> const ChemicalProperties&;

    /// Return the sensitivity of the equilibrium state.
    /// The sensitivity of the equilibrium state is defined as the rate of change of the
    /// molar amounts of the equilibrium species with respect to temperature `T`, pressure `P`,
    /// and molar amounts of equilibrium elements `be`.
    auto sensitivity() -> const EquilibriumSensitivity&;

    /// Compute the sensitivity of the species amounts with respect to temperature.
    auto dndT() -> VectorConstRef;

    /// Compute the sensitivity of the species amounts with respect to pressure.
    auto dndP() -> VectorConstRef;

    /// Compute the sensitivity of the species amounts with respect to element amounts.
    auto dndb() -> VectorConstRef;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
