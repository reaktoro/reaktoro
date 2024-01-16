// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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
#include <Reaktoro/Common/HashUtils.hpp>
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConditions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumPredictor.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/ODML/ClusterConnectivity.hpp>
#include <Reaktoro/ODML/PriorityQueue.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalProps;
class ChemicalState;
class ChemicalSystem;
class EquilibriumRestrictions;
class EquilibriumSpecs;
struct SmartEquilibriumOptions;
struct SmartEquilibriumResult;

/// Used for calculating chemical equilibrium states using an on-demand machine learning (ODML) strategy.
class SmartEquilibriumSolver
{
public:
    /// Construct an SmartEquilibriumSolver object with given chemical system.
    explicit SmartEquilibriumSolver(ChemicalSystem const& system);

    /// Construct an SmartEquilibriumSolver object with given chemical equilibrium specifications.
    explicit SmartEquilibriumSolver(EquilibriumSpecs const& specs);

    /// Construct a copy of an SmartEquilibriumSolver object.
    SmartEquilibriumSolver(SmartEquilibriumSolver const& other);

    /// Destroy this SmartEquilibriumSolver object.
    ~SmartEquilibriumSolver();

    /// Assign a copy of an SmartEquilibriumSolver object to this.
    auto operator=(SmartEquilibriumSolver other) -> SmartEquilibriumSolver&;

    //=================================================================================================================
    //
    // CHEMICAL EQUILIBRIUM METHODS
    //
    //=================================================================================================================

    /// Equilibrate a chemical state.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed equilibrium state (out)
    auto solve(ChemicalState& state) -> SmartEquilibriumResult;

    /// Equilibrate a chemical state respecting given reactivity restrictions.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed equilibrium state (out)
    /// @param restrictions The reactivity restrictions on the amounts of selected species
    auto solve(ChemicalState& state, EquilibriumRestrictions const& restrictions) -> SmartEquilibriumResult;

    /// Equilibrate a chemical state respecting given constraint conditions.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed equilibrium state (out)
    /// @param conditions The specified constraint conditions to be attained at chemical equilibrium
    auto solve(ChemicalState& state, EquilibriumConditions const& conditions) -> SmartEquilibriumResult;

    /// Equilibrate a chemical state respecting given constraint conditions and reactivity restrictions.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed equilibrium state (out)
    /// @param conditions The specified constraint conditions to be attained at chemical equilibrium
    /// @param restrictions The reactivity restrictions on the amounts of selected species
    auto solve(ChemicalState& state, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions) -> SmartEquilibriumResult;

    //=================================================================================================================
    //
    // CHEMICAL EQUILIBRIUM METHODS WITH SENSITIVITY CALCULATION
    //
    //=================================================================================================================

    /// Equilibrate a chemical state and compute sensitivity derivatives.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed equilibrium state (out)
    /// @param[out] sensitivity The sensitivity derivatives of the equilibrium state with respect to given input conditions
    auto solve(ChemicalState& state, EquilibriumSensitivity& sensitivity) -> SmartEquilibriumResult;

    /// Equilibrate a chemical state respecting given reactivity restrictions and compute sensitivity derivatives.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed equilibrium state (out)
    /// @param[out] sensitivity The sensitivity derivatives of the equilibrium state with respect to given input conditions
    /// @param restrictions The reactivity restrictions on the amounts of selected species
    auto solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, EquilibriumRestrictions const& restrictions) -> SmartEquilibriumResult;

    /// Equilibrate a chemical state respecting given constraint conditions and compute sensitivity derivatives.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed equilibrium state (out)
    /// @param[out] sensitivity The sensitivity derivatives of the equilibrium state with respect to given input conditions
    /// @param conditions The specified constraint conditions to be attained at chemical equilibrium
    auto solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, EquilibriumConditions const& conditions) -> SmartEquilibriumResult;

    /// Equilibrate a chemical state respecting given constraint conditions and reactivity restrictions and compute sensitivity derivatives.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed equilibrium state (out)
    /// @param[out] sensitivity The sensitivity derivatives of the equilibrium state with respect to given input conditions
    /// @param conditions The specified constraint conditions to be attained at chemical equilibrium
    /// @param restrictions The reactivity restrictions on the amounts of selected species
    auto solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions) -> SmartEquilibriumResult;

    //=================================================================================================================
    //
    // MISCELLANEOUS METHODS
    //
    //=================================================================================================================

    /// Set the options of the equilibrium solver.
    auto setOptions(SmartEquilibriumOptions const& options) -> void;

    /// The record of the knowledge database containing input, output, and derivatives data.
    struct Record
    {
        /// The fully calculated chemical equilibrium state.
        ChemicalState state;

        /// The conditions at which the chemical equilibrium state was calculated.
        EquilibriumConditions conditions;

        /// The sensitivity derivatives at the calculated chemical equilibrium state.
        EquilibriumSensitivity sensitivity;

        /// The predictor of chemical equilibrium states at given new conditions.
        EquilibriumPredictor predictor;
    };

    /// The cluster storing learned input-output data with same classification.
    struct Cluster
    {
        /// The indices of the primary species for this cluster.
        ArrayXl iprimary;

        /// The hash of the indices of the primary species for this cluster.
        Index label = 0;

        /// The records stored in this cluster with learning data.
        Deque<Record> records;

        /// The priority queue for the records based on their usage count.
        PriorityQueue priority;
    };

    /// The collection of clusters containing learned input-output data associated to a temperature-pressure grid cell.
    struct Cell
    {
        /// The clusters containing the learned input-output data points in a temperature-pressure grid cell.
        Deque<Cluster> clusters;

        /// The connectivity matrix of the clusters to determine how we move from one to another when searching.
        ClusterConnectivity connectivity;

        /// The priority queue for the clusters based on their usage counts.
        PriorityQueue priority;
    };

    /// The temperature-pressure grid cells containing learned input-output data.
    struct Grid
    {
        /// The hash table used to access a temperature-pressure grid cell containing learned computations.
        /// Note the use of `long` as number type for a temperature-pressure pairs. Temperatures and
        /// pressures are rounded to nearest checkpoints based on provided temperature/pressure step
        /// lengths for discretization.
        Map<Pair<long, long>, Cell> cells;
    };

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

} // namespace Reaktoro
