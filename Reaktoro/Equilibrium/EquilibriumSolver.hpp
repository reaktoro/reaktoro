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

#pragma once

// C++ includes
#include <memory>

// Reaktoro includes
#include <Reaktoro/Common/Matrix.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalState;
class ChemicalSystem;
class Partition;
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

    /// Find an initial feasible guess for an equilibrium problem
    /// @param state[in,out] The initial guess and the final state of the equilibrium approximation
    /// @param be The molar amounts of the elements in the equilibrium partition
    auto approximate(ChemicalState& state, double T, double P, const Vector& be) -> EquilibriumResult;

    /// Solve an equilibrium problem with given molar amounts of the elements in the equilibrium partition..
    /// @param state[in,out] The initial guess and the final state of the equilibrium calculation
    /// @param be The molar amounts of the elements in the equilibrium partition
    auto solve(ChemicalState& state, double T, double P, const Vector& be) -> EquilibriumResult;

    /// Solve an equilibrium problem with given molar amounts of the elements in the equilibrium partition..
    /// @param state[in,out] The initial guess and the final state of the equilibrium calculation
    /// @param T The temperature (in units of K)
    /// @param P The pressure (in units of Pa)
    /// @param be The molar amounts of the elements in the equilibrium partition
    auto solve(ChemicalState& state, double T, double P, const double* be) -> EquilibriumResult;

    /// Return the sensitivity of the equilibrium state.
    /// The sensitivity of the equilibrium state is defined as the rate of change of the
    /// molar amounts of the equilibrium species with respect to temperature `T`, pressure `P`,
    /// and molar amounts of equilibrium elements `be`.
    auto sensitivity() -> EquilibriumSensitivity;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
