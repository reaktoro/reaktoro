// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalProps;
class ChemicalState;
class ChemicalSystem;
class EquilibriumConditions;
class EquilibriumRestrictions;
class EquilibriumSensitivity;
class EquilibriumSpecs;
struct EquilibriumOptions;
struct EquilibriumResult;

/// A solver class for solving chemical equilibrium calculations.
class EquilibriumSolver
{
public:
    /// Construct an EquilibriumSolver object with given chemical system.
    explicit EquilibriumSolver(const ChemicalSystem& system);

    /// Construct an EquilibriumSolver object with given chemical equilibrium specifications.
    explicit EquilibriumSolver(const EquilibriumSpecs& specs);

    /// Construct a copy of an EquilibriumSolver object.
    EquilibriumSolver(const EquilibriumSolver& other);

    /// Destroy this EquilibriumSolver object.
    ~EquilibriumSolver();

    /// Assign a copy of an EquilibriumSolver object to this.
    auto operator=(EquilibriumSolver other) -> EquilibriumSolver&;

    /// Set the options of the equilibrium solver.
    auto setOptions(const EquilibriumOptions& options) -> void;

    /// Solve an equilibrium problem with given chemical state in disequilibrium.
    /// @param state[in,out] The initial guess and the final state of the equilibrium calculation
    auto solve(ChemicalState& state) -> EquilibriumResult;

    /// Solve an equilibrium problem with given chemical state in disequilibrium.
    /// @param state[in,out] The initial guess and the final state of the equilibrium calculation
    /// @param restrictions The restrictions on the reactivity amounts of the species
    auto solve(ChemicalState& state, const EquilibriumRestrictions& restrictions) -> EquilibriumResult;

    /// Solve an equilibrium problem with given chemical state in disequilibrium and equilibrium conditions.
    /// @param state[in,out] The initial guess and the final state of the equilibrium calculation
    /// @param conditions The conditions to be attained at chemical equilibrium
    auto solve(ChemicalState& state, const EquilibriumConditions& conditions) -> EquilibriumResult;

    /// Solve an equilibrium problem with given chemical state in disequilibrium and equilibrium conditions.
    /// @param state[in,out] The initial guess and the final state of the equilibrium calculation
    /// @param conditions The conditions to be attained at chemical equilibrium
    /// @param restrictions The restrictions on the reactivity amounts of the species
    auto solve(ChemicalState& state, const EquilibriumConditions& conditions, const EquilibriumRestrictions& restrictions) -> EquilibriumResult;

    /// Solve an equilibrium problem with given chemical state in disequilibrium.
    /// @param state[in,out] The initial guess and the final state of the equilibrium calculation
    /// @param sensitivity[out] The sensitivity of the equilibrium state with respect to its input conditions
    auto solve(ChemicalState& state, EquilibriumSensitivity& sensitivity) -> EquilibriumResult;

    /// Solve an equilibrium problem with given chemical state in disequilibrium.
    /// @param state[in,out] The initial guess and the final state of the equilibrium calculation
    /// @param sensitivity[out] The sensitivity of the equilibrium state with respect to its input conditions
    /// @param restrictions The restrictions on the reactivity amounts of the species
    auto solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, const EquilibriumRestrictions& restrictions) -> EquilibriumResult;

    /// Solve an equilibrium problem with given chemical state in disequilibrium and equilibrium conditions.
    /// @param state[in,out] The initial guess and the final state of the equilibrium calculation
    /// @param sensitivity[out] The sensitivity of the equilibrium state with respect to its input conditions
    /// @param conditions The conditions to be attained at chemical equilibrium
    auto solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, const EquilibriumConditions& conditions) -> EquilibriumResult;

    /// Solve an equilibrium problem with given chemical state in disequilibrium and equilibrium conditions.
    /// @param state[in,out] The initial guess and the final state of the equilibrium calculation
    /// @param sensitivity[out] The sensitivity of the equilibrium state with respect to its input conditions
    /// @param conditions The conditions to be attained at chemical equilibrium
    /// @param restrictions The restrictions on the reactivity amounts of the species
    auto solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, const EquilibriumConditions& conditions, const EquilibriumRestrictions& restrictions) -> EquilibriumResult;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

} // namespace Reaktoro
