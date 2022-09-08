// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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
#include <Reaktoro/Common/Matrix.hpp>

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

/// Used for calculating chemical equilibrium states.
class EquilibriumSolver
{
public:
    /// Construct an EquilibriumSolver object with given chemical system.
    explicit EquilibriumSolver(ChemicalSystem const& system);

    /// Construct an EquilibriumSolver object with given chemical equilibrium specifications.
    explicit EquilibriumSolver(EquilibriumSpecs const& specs);

    /// Construct a copy of an EquilibriumSolver object.
    EquilibriumSolver(EquilibriumSolver const& other);

    /// Destroy this EquilibriumSolver object.
    ~EquilibriumSolver();

    /// Assign a copy of an EquilibriumSolver object to this.
    auto operator=(EquilibriumSolver other) -> EquilibriumSolver&;

    //=================================================================================================================
    //
    // CHEMICAL EQUILIBRIUM METHODS
    //
    //=================================================================================================================

    /// Equilibrate a chemical state.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed equilibrium state (out)
    auto solve(ChemicalState& state) -> EquilibriumResult;

    /// Equilibrate a chemical state respecting given reactivity restrictions.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed equilibrium state (out)
    /// @param restrictions The reactivity restrictions on the amounts of selected species
    auto solve(ChemicalState& state, EquilibriumRestrictions const& restrictions) -> EquilibriumResult;

    /// Equilibrate a chemical state respecting given constraint conditions.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed equilibrium state (out)
    /// @param conditions The specified constraint conditions to be attained at chemical equilibrium
    auto solve(ChemicalState& state, EquilibriumConditions const& conditions) -> EquilibriumResult;

    /// Equilibrate a chemical state respecting given constraint conditions and reactivity restrictions.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed equilibrium state (out)
    /// @param conditions The specified constraint conditions to be attained at chemical equilibrium
    /// @param restrictions The reactivity restrictions on the amounts of selected species
    auto solve(ChemicalState& state, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions) -> EquilibriumResult;

    //=================================================================================================================
    //
    // CHEMICAL EQUILIBRIUM METHODS WITH SENSITIVITY CALCULATION
    //
    //=================================================================================================================

    /// Equilibrate a chemical state and compute sensitivity derivatives.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed equilibrium state (out)
    /// @param[out] sensitivity The sensitivity derivatives of the equilibrium state with respect to given input conditions
    auto solve(ChemicalState& state, EquilibriumSensitivity& sensitivity) -> EquilibriumResult;

    /// Equilibrate a chemical state respecting given reactivity restrictions and compute sensitivity derivatives.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed equilibrium state (out)
    /// @param[out] sensitivity The sensitivity derivatives of the equilibrium state with respect to given input conditions
    /// @param restrictions The reactivity restrictions on the amounts of selected species
    auto solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, EquilibriumRestrictions const& restrictions) -> EquilibriumResult;

    /// Equilibrate a chemical state respecting given constraint conditions and compute sensitivity derivatives.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed equilibrium state (out)
    /// @param[out] sensitivity The sensitivity derivatives of the equilibrium state with respect to given input conditions
    /// @param conditions The specified constraint conditions to be attained at chemical equilibrium
    auto solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, EquilibriumConditions const& conditions) -> EquilibriumResult;

    /// Equilibrate a chemical state respecting given constraint conditions and reactivity restrictions and compute sensitivity derivatives.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed equilibrium state (out)
    /// @param[out] sensitivity The sensitivity derivatives of the equilibrium state with respect to given input conditions
    /// @param conditions The specified constraint conditions to be attained at chemical equilibrium
    /// @param restrictions The reactivity restrictions on the amounts of selected species
    auto solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions) -> EquilibriumResult;

    //=================================================================================================================
    //
    // MISCELLANEOUS METHODS
    //
    //=================================================================================================================

    /// Set the options of the equilibrium solver.
    auto setOptions(EquilibriumOptions const& options) -> void;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

} // namespace Reaktoro
