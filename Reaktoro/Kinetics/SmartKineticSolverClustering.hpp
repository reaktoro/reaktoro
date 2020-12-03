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
#include <string>

// Reaktoro includes
#include <Reaktoro/Math/Matrix.hpp>
#include <Reaktoro/Common/ChemicalVector.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Profiling.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumSolver.hpp>
#include <Reaktoro/Kinetics/SmartKineticOptions.hpp>
#include <Reaktoro/Kinetics/SmartKineticResult.hpp>
#include <Reaktoro/ODML/ClusterConnectivity.hpp> // class needed to maintain connectivity of the clusters
#include <Reaktoro/ODML/PriorityQueue.hpp> // class needed to maintain priority queues of clusters and records inside them

#include <Reaktoro/Kinetics/SmartKineticSolverBase.hpp>

namespace Reaktoro {

/// A class that represents a solver for chemical kinetics problems.
/// @see KineticProblem
class SmartKineticSolverClustering : public SmartKineticSolverBase
{
public:
    /// Construct a SmartKineticSolverClustering instance.
    explicit SmartKineticSolverClustering(const ReactionSystem& reactions, const Partition& partition);

    /// Destroy the SmartKineticSolverClustering instance.
    virtual ~SmartKineticSolverClustering();

    /// Learn how to perform a full equilibrium calculation (with tracking)
    auto learn(ChemicalState& state, double& t, double dt) -> void;

    /// Estimate the equilibrium state using sensitivity derivatives (profiling the expenses)
    auto estimate(ChemicalState& state, double& t, double dt) -> void;

    /// Output clusters created during the ODML algorithm.
    auto outputInfo() const -> void;

private:

    /// A structure used to store the record of the knowledge database containing input, output, and derivatives data.
    struct KineticRecord
    {
        /// The temperature of the equilibrium state (in units of K).
        double T;

        /// The pressure of the equilibrium state (in units of Pa).
        double P;

        /// The amounts of elements in the equilibrium state (in units of mol).
        Vector be;

        /// The calculated equilibrium state at `T`, `P`, `be`.
        ChemicalState state;

        /// The chemical properties at the calculated equilibrium state.
        ChemicalProperties properties;

        /// The sensitivity derivatives at the calculated equilibrium state.
        EquilibriumSensitivity sensitivity;

        /// The matrix used to compute relative change of chemical potentials due to change in `be`.
        Matrix Mbe;

        /// Calculated kinetic state, containing initial and final value of [be, nk], initial and final time, and sensitivities
        ODEState ode_state;

        // Chemical kinetic rates
        ChemicalVector rates;
    };

    /// The cluster storing learned input-output data with same classification.
    struct Cluster
    {
        /// The indices of the primary species for this cluster.
        VectorXi iprimary;

        /// The hash of the indices of the primary species for this cluster.
        std::size_t label = 0;

        /// The records stored in this cluster with learning data.
        std::deque<KineticRecord> records;

        /// The priority queue for the records based on their usage count.
        PriorityQueue priority;
    };

    /// The database containing the learned input-output data.
    struct Database
    {
        /// The clusters containing the learned input-output data points.
        std::deque<Cluster> clusters;

        /// The connectivity matrix of the clusters to determine how we move from one to another when searching.
        ClusterConnectivity connectivity;

        /// The priority queue for the clusters based on their usage counts.
        PriorityQueue priority;
    };

    /// The database with learned input-output data points.
    Database database;

    /// The storage for matrix du/db = du/dn * dn/db
    Matrix dudb;

    /// The storage for vector u(iprimary)
    Vector up;

    /// The storage for matrix Mbe = inv(u(iprimary)) * du(iprimary)/db.
    Matrix Mbe;

};

/// A class that represents a solver for chemical kinetics problems.
/// @see KineticProblem
class SmartKineticSolverClusteringExtended : public SmartKineticSolverBase
{
public:
    /// Construct a SmartKineticSolverClustering instance.
    explicit SmartKineticSolverClusteringExtended(const ReactionSystem& reactions, const Partition& partition);

    /// Destroy the SmartKineticSolverClustering instance.
    virtual ~SmartKineticSolverClusteringExtended();

    /// Learn how to perform a full equilibrium calculation (with tracking)
    auto learn(ChemicalState& state, double& t, double dt) -> void;

    /// Estimate the equilibrium state using sensitivity derivatives (profiling the expenses)
    auto estimate(ChemicalState& state, double& t, double dt) -> void;

    /// Output clusters created during the ODML algorithm.
    auto outputInfo() const -> void;

private:

    /// A structure used to store the record of the knowledge database containing input, output, and derivatives data.
    struct KineticRecord
    {
        /// The temperature of the equilibrium state (in units of K).
        double T;

        /// The pressure of the equilibrium state (in units of Pa).
        double P;

        /// The amounts of elements in the equilibrium state (in units of mol).
        Vector be;

        /// The calculated equilibrium state at `T`, `P`, `be`.
        ChemicalState state;

        /// The chemical properties at the calculated equilibrium state.
        ChemicalProperties properties;

        /// The sensitivity derivatives at the calculated equilibrium state.
        EquilibriumSensitivity sensitivity;

        /// The matrix used to compute relative change of chemical potentials due to change in `be`.
        Matrix Mbe;

        /// Calculated kinetic state, containing initial and final value of [be, nk], initial and final time, and sensitivities
        ODEState ode_state;

        // Chemical kinetic rates
        ChemicalVector rates;
    };

    /// The cluster (extended to kinetics) storing learned input-output data with same classification.
    struct ClusterExtended
    {
        /// The indices of the primary species for this cluster.
        VectorXi iprimary;

        /// The indices of the primary species and kinetics species for this cluster.
        VectorXi iprimary_ikin;

        /// The hash of the indices of the primary species for this cluster.
        std::size_t label = 0;

        /// The hash of the indices of the primary and kinetics species for this cluster.
        std::size_t label_extended = 0;

        /// The records stored in this cluster with learning data.
        std::deque<KineticRecord> records;

        /// The priority queue for the records based on their usage count.
        PriorityQueue priority;
    };

    /// The database containing the learned input-output data.
    struct DatabaseExtended
    {
        /// The clusters containing the learned input-output data points.
        std::deque<ClusterExtended> clusters;

        /// The connectivity matrix of the clusters to determine how we move from one to another when searching.
        ClusterConnectivity connectivity;

        /// The priority queue for the clusters based on their usage counts.
        PriorityQueue priority;
    };

    /// The database with learned input-output data points.
    DatabaseExtended database_extended;

    /// The storage for matrix du/db = du/dn * dn/db
    Matrix dudb;

    /// The storage for vector u(iprimary)
    Vector up;

    /// The storage for matrix Mbe = inv(u(iprimary)) * du(iprimary)/db.
    Matrix Mbe;

};

} // namespace Reaktoro
