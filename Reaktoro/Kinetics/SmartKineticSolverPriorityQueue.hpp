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
#include <deque>

// Reaktoro includes
#include <Reaktoro/Math/Matrix.hpp>
#include <Reaktoro/Common/ChemicalVector.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Kinetics/SmartKineticSolverBase.hpp>

namespace Reaktoro {

/// A class that represents a solver for chemical kinetics problems.
/// @see KineticProblem
class SmartKineticSolverPriorityQueue : public SmartKineticSolverBase
{
public:
    /// Construct a SmartKineticSolverPriorityQueue instance.
    explicit SmartKineticSolverPriorityQueue(const ReactionSystem& reactions, const Partition& partition);

    /// Destroy the SmartKineticSolverPriorityQueue instance.
    virtual ~SmartKineticSolverPriorityQueue();

    /// Learn how to perform a full equilibrium calculation (with tracking)
    auto learn(ChemicalState& state, double& t, double dt) -> void;

    /// Estimate the equilibrium state using sensitivity derivatives (profiling the expenses)
    auto estimate(ChemicalState& state, double& t, double dt) -> void;

    /// Output clusters created during the ODML algorithm.
    auto outputInfo() const -> void;

private:

    /// A structure used to store the node of tree for smart kinetic calculations.
    struct KineticRecord{

        /// The calculated equilibrium state at `T`, `P`, `be`.
        ChemicalState chemical_state;

        /// The chemical properties at the calculated equilibrium state.
        ChemicalProperties properties;

        /// The sensitivity derivatives at the calculated equilibrium state.
        EquilibriumSensitivity sensitivity;

        /// The matrix used to compute relative change of chemical potentials due to change in `be`.
        Matrix Mb;

        // Indices of the major species
        VectorXi imajor;

        /// Calculated kinetic state, containing initial and final value of [be, nk], initial and final time, and sensitivities
        ODEState ode_state;

        // Chemical kinetic rates
        ChemicalVector rates;
    };

    /// The tree used to save the calculated equilibrium states and respective sensitivities
    std::deque<KineticRecord> database;

    // The queue with priority indices, in which equilibrium states of the `tree` must be considered
    std::deque<Index> kinetics_priority;

    // The queue with the rank of each equilibrium states in the `tree`
    std::deque<Index> kinetics_ranking;

    /// The chemical potentials at the calculated equilibrium state
    ChemicalVector u;

    /// The storage for matrix du/db = du/dn * dn/db
    Matrix dudb;

    /// The storage for vector u(imajor)
    Vector um;

    /// The storage for matrix Mbe = inv(u(imajor)) * du(imajor)/db.
    Matrix Mbe;

};

/// A class that represents a solver for chemical kinetics problems.
/// @see KineticProblem
class SmartKineticSolverPriorityQueuePrimary : public SmartKineticSolverBase
{
public:
    /// Construct a SmartKineticSolverPriorityQueue instance.
    explicit SmartKineticSolverPriorityQueuePrimary(const ReactionSystem& reactions, const Partition& partition);

    /// Destroy the SmartKineticSolverPriorityQueue instance.
    virtual ~SmartKineticSolverPriorityQueuePrimary();

    /// Learn how to perform a full equilibrium calculation (with tracking)
    auto learn(ChemicalState& state, double& t, double dt) -> void;

    /// Estimate the equilibrium state using sensitivity derivatives (profiling the expenses)
    auto estimate(ChemicalState& state, double& t, double dt) -> void;

    /// Output clusters created during the ODML algorithm.
    auto outputInfo() const -> void;

private:

    /// A structure used to store the node of tree for smart kinetic calculations.
    struct KineticRecord{

        /// The calculated equilibrium state at `T`, `P`, `be`.
        ChemicalState chemical_state;

        /// The chemical properties at the calculated equilibrium state.
        ChemicalProperties properties;

        /// The sensitivity derivatives at the calculated equilibrium state.
        EquilibriumSensitivity sensitivity;

        /// The matrix used to compute relative change of chemical potentials due to change in `be`.
        Matrix Mb;

        // Indices of the major species
        VectorXi imajor;

        /// Calculated kinetic state, containing initial and final value of [be, nk], initial and final time, and sensitivities
        ODEState ode_state;

        // Chemical kinetic rates
        ChemicalVector rates;
    };

    /// The tree used to save the calculated equilibrium states and respective sensitivities
    std::deque<KineticRecord> database;

    // The queue with priority indices, in which equilibrium states of the `tree` must be considered
    std::deque<Index> kinetics_priority;

    // The queue with the rank of each equilibrium states in the `tree`
    std::deque<Index> kinetics_ranking;

    /// The chemical potentials at the calculated equilibrium state
    ChemicalVector u;

    /// The storage for matrix du/db = du/dn * dn/db
    Matrix dudb;

    /// The storage for vector u(iprimary)
    Vector up;

    /// The storage for matrix Mbe = inv(u(iprimary)) * du(iprimary)/db.
    Matrix Mbe;

};

} // namespace Reaktoro
