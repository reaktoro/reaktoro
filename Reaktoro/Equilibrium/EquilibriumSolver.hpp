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
class EquilibriumProblem;
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
    explicit EquilibriumSolver(const ChemicalSystem& system);

    /// Construct an EquilibriumSolver object with given chemical equilibrium specifications.
    explicit EquilibriumSolver(const EquilibriumSpecs& specs);

    /// Construct a copy of an EquilibriumSolver object.
    EquilibriumSolver(const EquilibriumSolver& other);

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
    auto solve(ChemicalState& state, const EquilibriumRestrictions& restrictions) -> EquilibriumResult;

    /// Equilibrate a chemical state respecting given constraint conditions.
    /// Solve an equilibrium problem with given chemical state and equilibrium conditions.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed equilibrium state (out)
    /// @param conditions The specified constraint conditions to be attained at chemical equilibrium
    auto solve(ChemicalState& state, const EquilibriumConditions& conditions) -> EquilibriumResult;

    /// Equilibrate a chemical state respecting given constraint conditions and reactivity restrictions.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed equilibrium state (out)
    /// @param conditions The specified constraint conditions to be attained at chemical equilibrium
    /// @param restrictions The reactivity restrictions on the amounts of selected species
    auto solve(ChemicalState& state, const EquilibriumConditions& conditions, const EquilibriumRestrictions& restrictions) -> EquilibriumResult;

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
    auto solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, const EquilibriumRestrictions& restrictions) -> EquilibriumResult;

    /// Equilibrate a chemical state respecting given constraint conditions and compute sensitivity derivatives.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed equilibrium state (out)
    /// @param[out] sensitivity The sensitivity derivatives of the equilibrium state with respect to given input conditions
    /// @param conditions The specified constraint conditions to be attained at chemical equilibrium
    auto solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, const EquilibriumConditions& conditions) -> EquilibriumResult;

    /// Equilibrate a chemical state respecting given constraint conditions and reactivity restrictions and compute sensitivity derivatives.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed equilibrium state (out)
    /// @param[out] sensitivity The sensitivity derivatives of the equilibrium state with respect to given input conditions
    /// @param conditions The specified constraint conditions to be attained at chemical equilibrium
    /// @param restrictions The reactivity restrictions on the amounts of selected species
    auto solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, const EquilibriumConditions& conditions, const EquilibriumRestrictions& restrictions) -> EquilibriumResult;

    //=================================================================================================================
    //
    // CHEMICAL EQUILIBRIUM METHODS WITH GIVEN AMOUNTS OF CONSERVATIVE COMPONENTS
    //
    //=================================================================================================================

    /// Equilibrate a chemical state.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed equilibrium state (out)
    /// @param b0 The amounts of the conservative components in the chemical equilibrium problem
    auto solve(ChemicalState& state, ArrayXdConstRef b0) -> EquilibriumResult;

    /// Equilibrate a chemical state respecting given reactivity restrictions.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed equilibrium state (out)
    /// @param restrictions The reactivity restrictions on the amounts of selected species
    /// @param b0 The amounts of the conservative components in the chemical equilibrium problem
    auto solve(ChemicalState& state, const EquilibriumRestrictions& restrictions, ArrayXdConstRef b0) -> EquilibriumResult;

    /// Equilibrate a chemical state respecting given constraint conditions.
    /// Solve an equilibrium problem with given chemical state and equilibrium conditions.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed equilibrium state (out)
    /// @param conditions The specified constraint conditions to be attained at chemical equilibrium
    /// @param b0 The amounts of the conservative components in the chemical equilibrium problem
    auto solve(ChemicalState& state, const EquilibriumConditions& conditions, ArrayXdConstRef b0) -> EquilibriumResult;

    /// Equilibrate a chemical state respecting given constraint conditions and reactivity restrictions.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed equilibrium state (out)
    /// @param conditions The specified constraint conditions to be attained at chemical equilibrium
    /// @param restrictions The reactivity restrictions on the amounts of selected species
    /// @param b0 The amounts of the conservative components in the chemical equilibrium problem
    auto solve(ChemicalState& state, const EquilibriumConditions& conditions, const EquilibriumRestrictions& restrictions, ArrayXdConstRef b0) -> EquilibriumResult;

    //=================================================================================================================
    //
    // CHEMICAL EQUILIBRIUM METHODS WITH GIVEN AMOUNTS OF CONSERVATIVE COMPONENTS AND SENSITIVITY CALCULATION
    //
    //=================================================================================================================

    /// Equilibrate a chemical state and compute sensitivity derivatives.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed equilibrium state (out)
    /// @param[out] sensitivity The sensitivity derivatives of the equilibrium state with respect to given input conditions
    /// @param b0 The amounts of the conservative components in the chemical equilibrium problem
    auto solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, ArrayXdConstRef b0) -> EquilibriumResult;

    /// Equilibrate a chemical state respecting given reactivity restrictions and compute sensitivity derivatives.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed equilibrium state (out)
    /// @param[out] sensitivity The sensitivity derivatives of the equilibrium state with respect to given input conditions
    /// @param restrictions The reactivity restrictions on the amounts of selected species
    /// @param b0 The amounts of the conservative components in the chemical equilibrium problem
    auto solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, const EquilibriumRestrictions& restrictions, ArrayXdConstRef b0) -> EquilibriumResult;

    /// Equilibrate a chemical state respecting given constraint conditions and compute sensitivity derivatives.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed equilibrium state (out)
    /// @param[out] sensitivity The sensitivity derivatives of the equilibrium state with respect to given input conditions
    /// @param conditions The specified constraint conditions to be attained at chemical equilibrium
    /// @param b0 The amounts of the conservative components in the chemical equilibrium problem
    auto solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, const EquilibriumConditions& conditions, ArrayXdConstRef b0) -> EquilibriumResult;

    /// Equilibrate a chemical state respecting given constraint conditions and reactivity restrictions and compute sensitivity derivatives.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed equilibrium state (out)
    /// @param[out] sensitivity The sensitivity derivatives of the equilibrium state with respect to given input conditions
    /// @param conditions The specified constraint conditions to be attained at chemical equilibrium
    /// @param restrictions The reactivity restrictions on the amounts of selected species
    /// @param b0 The amounts of the conservative components in the chemical equilibrium problem
    auto solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, const EquilibriumConditions& conditions, const EquilibriumRestrictions& restrictions, ArrayXdConstRef b0) -> EquilibriumResult;

    //=================================================================================================================
    //
    // MISCELLANEOUS METHODS
    //
    //=================================================================================================================

    /// Set the options of the equilibrium solver.
    auto setOptions(const EquilibriumOptions& options) -> void;

    /// Return the conservative matrix of the chemical equilibrium problem.
    /// This conservative matrix of the chemical equilibrium problem is a
    /// matrix whose upper rows contains the formula matrix of the species with
    /// respect to elements and electric charge, and the lower rows contains
    /// the coefficients of the reactivity constraints (e.g., the
    /// stoichiometric matrix of the restricted reactions in the equilibrium
    /// computation). This matrix is used to compute the amounts of the
    /// conservative components in the chemical equilibrium problem, which are
    /// elements, electric charge, and the extend of the restricted reactions.
    auto conservativeMatrix() const -> MatrixXdConstRef;

    /// Return the amounts of the conservative components (elements, charge, extent of restricted reactions) in the chemical equilibrium problem.
    auto componentAmounts(ChemicalState const& state) const -> ArrayXr;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

} // namespace Reaktoro
