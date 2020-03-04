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

#include "SmartEquilibriumSolver.hpp"

// C++ includes
#include <algorithm>
#include <deque>
#include <numeric>
#include <set>
#include <tuple>

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Profiling.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumResult.hpp>
#include <Reaktoro/ODML/ClusterConnectivity.hpp>
#include <Reaktoro/ODML/PriorityQueue.hpp>
#include <Reaktoro/Optimization/Canonicalizer.hpp>

namespace Reaktoro {

struct SmartEquilibriumSolver::Impl
{
    /// The chemical system instance
    ChemicalSystem system;

    /// The partition of the chemical system
    Partition partition;

    /// The canonicalizer used to determine primary and secondary species
    Canonicalizer canonicalizer;

    /// The options for the smart equilibrium calculation
    SmartEquilibriumOptions options;

    /// The result of the last smart equilibrium calculation.
    SmartEquilibriumResult result;

    /// The solver for the equilibrium calculations
    EquilibriumSolver solver;

    /// The chemical properties of the chemical system
    ChemicalProperties properties;

    /// The amounts of the species in the chemical system
    Vector n;

    /// The amounts of the elements in the equilibrium partition
    Vector be;

    /// The storage for matrix du/db = du/dn * dn/db
    Matrix dudb;

    /// The storage for vector u(iprimary)
    Vector up;

    /// The storage for matrix Mb = inv(u(iprimary)) * du(iprimary)/db.
    Matrix Mb;

    /// The record of the knowledge database containing input, output and derivatives data.
    struct Record
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
        Matrix Mb;
    };

    /// The cluster storing learned input-output data with same classification.
    struct Cluster
    {
        /// The indices of the primary species for this cluster.
        VectorXi iprimary;

        /// The records stored in this cluster with learning data.
        std::deque<Record> records;

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

    /// Construct a default SmartEquilibriumSolver::Impl instance.
    Impl()
    {}

    /// Construct an SmartEquilibriumSolver::Impl instance with given partition.
    Impl(const Partition& partition_)
    : partition(partition_), system(partition_.system()), solver(partition_)
    {
        // Initialize the canonicalizer with the formula matrix A
        canonicalizer.compute(system.formulaMatrix());
    }

    /// Set the options for the equilibrium calculation.
    auto setOptions(const SmartEquilibriumOptions& options_) -> void
    {
        options = options_;

        // Tweak the options for the Gibbs energy minimization during learning operations.
        options.learning.hessian = GibbsHessian::Exact; // ensure the use of an exact Hessian of the Gibbs energy function
        options.learning.optimum.tolerance = 1e-10; // ensure the use of a stricter residual tolerance for the Gibbs energy minimization

        solver.setOptions(options.learning);
    }

    /// Set the partition of the chemical system.
    auto setPartition(const Partition& partition_) -> void
    {
        partition = partition_;

        solver.setPartition(partition_);
    }

    /// Learn how to perform a full equilibrium calculation (with tracking)
    auto learn(ChemicalState& state, double T, double P, VectorConstRef be) -> void
    {
        // Calculate the equilibrium state using conventional Gibbs energy minimization approach
        timeit( solver.solve(state, T, P, be),
            result.timing.learning_gibbs_energy_minimization= );

        // Store the result of the Gibbs energy minimization calculation performed during learning
        result.learning.gibbs_energy_minimization = solver.result();

        // Update `properties` with the chemical properties of the calculated equilibrum state
        properties = solver.properties();

        // Get a reference to the sensitivity derivatives of the calculated equilibrium state
        const auto& sensitivity = solver.sensitivity();

        // The species amounts at the calculated equilibrium state
        const auto& n = state.speciesAmounts();

        // Update the canonical form of formula matrix A so that we can identify primary species
        canonicalizer.updateWithPriorityWeights(n);

        // Store the indices of primary and secondary species in state
        state.equilibrium().setIndicesEquilibriumSpecies(canonicalizer.Q(), canonicalizer.numBasicVariables());

        // The indices of the primary species at the calculated equilibrium state
        const auto& iprimary = state.equilibrium().indicesPrimarySpecies();

        // The chemical potentials at the calculated equilibrium state
        const auto u = properties.chemicalPotentials();

        // Auxiliary references to the derivatives dn/db and du/dn
        const auto& dndb = solver.sensitivity().dndb;
        const auto& dudn = u.ddn;

        // Compute the matrix du/db = du/dn * dn/db
        dudb = dudn * dndb;

        // The vector u(iprimary) with chemical potentials of primary species
        up.noalias() = rows(u.val, iprimary);

        // The matrix du(iprimary)/db with derivatives of chemical potentials of primary species
        const auto dupdb = rows(dudb, iprimary);

        // Compute matrix Mb = 1/up * dup/db
        Mb.noalias() = diag(inv(up)) * dupdb;

        // Find the cluster
        auto iter = std::find_if(database.clusters.begin(), database.clusters.end(),
            [&](const Cluster& cluster) { return cluster.iprimary == iprimary; });

        // If cluster found, store the new record in it, otherwise, create a new cluster
        if(iter < database.clusters.end())
        {
            auto& cluster = *iter;
            cluster.records.push_back({T, P, be, state, properties, sensitivity, Mb});
            cluster.priority.extend();
        }
        else
        {
            // Create a new cluster
            Cluster cluster;
            cluster.iprimary = iprimary;
            cluster.records.push_back({T, P, be, state, properties, sensitivity, Mb});
            cluster.priority.extend();

            // Append the new cluster in the database
            database.clusters.push_back(cluster);
            database.connectivity.extend();
            database.priority.extend();
        }

        // timeit( tree.push_back({be, be/sum(be), state, properties, solver.sensitivity(), Mb, iprimary}),
        //     result.timing.learning_storage= );
    }

    /// Estimate the equilibrium state using sensitivity derivatives (profiling the expences)
    auto estimate(ChemicalState& state, double T, double P, VectorConstRef be) -> void
    {
        result.estimate.accepted = false;

        // Skip estimation if no cluster exists yet
        if(database.clusters.empty())
            return;

        // Auxiliary relative and absolute tolerance parameters
        const auto reltol = options.reltol;
        const auto abstol = options.abstol;

        Vector dbe;

        // Check if an entry in the database pass the error test.
        // It returns (`success`, `error`, `ispecies`), where
        //   - `success` is true if error test succeeds, false otherwise.
        //   - `error` is the first error violating the tolerance
        //   - `ispecies` is the index of the species that fails the error test
        auto pass_error_test = [&](const Record& record) -> std::tuple<bool, double, Index>
        {
            using std::abs;
            using std::max;
            const auto& be0 = record.be;
            const auto& Mb0 = record.Mb;

            dbe.noalias() = be - be0;

            double error = 0.0;
            for(auto i = 0; i < Mb.rows(); ++i) {
                error = max(error, abs(Mb0.row(i) * dbe));
                if(error >= reltol)
                    return { false, error, i };
            }

            return { true, error, -1 };
        };

        // The current set of primary species in the chemical state
        const auto& iprimary = state.equilibrium().indicesPrimarySpecies();

        // The function that identifies the index of starting cluster
        auto index_starting_cluster = [&]() -> Index
        {
            // If no primary species, then return number of clusters to trigger use of total usage counts of clusters
            if(iprimary.size() == 0)
                return database.clusters.size();

            // Find the index of the cluster with same set of primary species (search those with highest count first)
            for(auto icluster : database.priority.order())
                if(database.clusters[icluster].iprimary == iprimary)
                    return icluster;

            return database.clusters.size();
        };

        // The index of the starting cluster
        const auto icluster = index_starting_cluster();

        // The ordering of the clusters to look for (starting with icluster)
        const auto& clusters_ordering = database.connectivity.order(icluster);

        // Iterate over all clusters starting with icluster
        for(auto jcluster : clusters_ordering)
        {
            const auto& records = database.clusters[jcluster].records;
            const auto& records_ordering = database.clusters[jcluster].priority.order();

            // Iterate over all records in current cluster
            for(auto irecord : records_ordering)
            {
                const auto& record = records[irecord];

                const auto [success, error, iprimaryspecies] = pass_error_test(record);

                if(success)
                {
                    const auto& be0 = record.be;
                    const auto& n0 = record.state.speciesAmounts();
                    const auto& y0 = record.state.equilibrium().elementChemicalPotentials();
                    const auto& z0 = record.state.equilibrium().speciesStabilities();
                    const auto& ies0 = record.state.equilibrium().indicesEquilibriumSpecies();
                    const auto& kp0 = record.state.equilibrium().numPrimarySpecies();
                    const auto& dndb0 = record.sensitivity.dndb;

                    n.noalias() = n0 + dndb0 * (be - be0);

                    const double nmin = min(n);  // TODO: this should be amounts of equilibrium species only!
                    const double nsum = sum(n);  // TODO: this should be amounts of equilibrium species only!

                    const auto eps_n = options.amount_fraction_cutoff * nsum;

                    // Check if all projected species amounts are not negative beyond tolerance
                    if(nmin < -eps_n)
                        continue;

                    // Increment priority of the current record (irecord) in the current cluster (jcluster)
                    database.clusters[jcluster].priority.increment(irecord);

                    // Increment priority of the current cluster (jcluster) with respect to starting cluster (icluster)
                    database.connectivity.increment(icluster, jcluster);

                    // Increment priority of the current cluster (jcluster)
                    database.priority.increment(jcluster);

                    // Set the chemical state result with estimated amounts
                    state.setSpeciesAmounts(n);
                    state.equilibrium().setElementChemicalPotentials(y0);
                    state.equilibrium().setSpeciesStabilities(z0);
                    state.equilibrium().setIndicesEquilibriumSpecies(ies0, kp0);

                    // Update the chemical properties of the system
                    properties = record.properties;  // TODO: We need to estimate properties = properties0 + variation : THIS IS A TEMPORARY SOLUTION!!!
                    result.estimate.accepted = true;
                    return;
                }
            }
        }

        result.estimate.accepted = false;
        return;
    }

    /// Solve the equilibrium problem with given initial state
    auto solve(ChemicalState& state) -> SmartEquilibriumResult
    {
        const auto& iee = partition.indicesEquilibriumElements();
        const auto& ies = partition.indicesEquilibriumSpecies();
        const auto T = state.temperature();
        const auto P = state.pressure();
        be = state.elementAmountsInSpecies(ies)(iee);
        return solve(state, T, P, be);
    }

    /// Solve the equilibrium problem with given problem definition
    auto solve(ChemicalState& state, const EquilibriumProblem& problem) -> SmartEquilibriumResult
    {
        setPartition(problem.partition());
        const auto T = problem.temperature();
        const auto P = problem.pressure();
        const auto& iee = partition.indicesEquilibriumElements(); // This statement needs to be after setPartition
        be = problem.elementAmounts()(iee);
        return solve(state, T, P, be);
    }

    auto solve(ChemicalState& state, double T, double P, VectorConstRef be) -> SmartEquilibriumResult
    {
        tic(0);

        // Absolutely ensure an exact Hessian of the Gibbs energy function is used in the calculations
        setOptions(options);

        // Reset the result of the last smart equilibrium calculation
        result = {};

        // Perform a smart estimate of the chemical state
        timeit( estimate(state, T, P, be),
            result.timing.estimate= );

        // Perform a learning step if the smart prediction is not sactisfatory
        if(!result.estimate.accepted)
            timeit( learn(state, T, P, be), result.timing.learn= );

        toc(0, result.timing.solve);

        return result;
    }

    auto outputClusterInfo() const -> void
    {
        for(auto i : database.priority.order())
        {
            std::cout << "CLUSTER #" << i << std::endl;
            std::cout << "  RANK OF CLUSTER: " << database.priority.priorities()[i] << std::endl;
            std::cout << "  PRIMARY SPECIES: ";
            for(auto j : database.clusters[i].iprimary)
                std::cout << system.species(j).name() << " ";
            std::cout << std::endl;
            std::cout << "  NUMBER OF RECORDS: " << database.clusters[i].records.size() << std::endl;
            std::cout << "  RANK OF RECORDS: ";
            for(auto j : database.clusters[i].priority.order())
                std::cout << database.clusters[i].priority.priorities()[j] << " ";
            std::cout << std::endl;
            std::cout << "  NEXT CLUSTER: " << database.connectivity.order(i)[1] << std::endl;
            std::cout << "  NEXT CLUSTER PRIMARY SPECIES: ";
            for(auto j : database.clusters[database.connectivity.order(i)[1]].iprimary)
                std::cout << system.species(j).name() << " ";
            std::cout << std::endl;
            std::cout << std::endl;
        }
    }
};

SmartEquilibriumSolver::SmartEquilibriumSolver()
: pimpl(new Impl())
{}

SmartEquilibriumSolver::SmartEquilibriumSolver(const ChemicalSystem& system)
: pimpl(new Impl(Partition(system)))
{}

SmartEquilibriumSolver::SmartEquilibriumSolver(const Partition& partition)
: pimpl(new Impl(partition))
{}

SmartEquilibriumSolver::SmartEquilibriumSolver(const SmartEquilibriumSolver& other)
: pimpl(new Impl(*other.pimpl))
{}

auto SmartEquilibriumSolver::operator=(SmartEquilibriumSolver other) -> SmartEquilibriumSolver&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

SmartEquilibriumSolver::~SmartEquilibriumSolver()
{
    // TODO Remove this from the desctructor.
    // TODO Remove this from the desctructor.
    // TODO Remove this from the desctructor.
    // TODO Remove this from the desctructor.
    // TODO Remove this from the desctructor.
    pimpl->outputClusterInfo();
}

auto SmartEquilibriumSolver::setOptions(const SmartEquilibriumOptions& options) -> void
{
    pimpl->setOptions(options);
}

auto SmartEquilibriumSolver::setPartition(const Partition& partition) -> void
{
    pimpl->setPartition(partition);
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, double T, double P, VectorConstRef be) -> SmartEquilibriumResult
{
    return pimpl->solve(state, T, P, be);
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, const EquilibriumProblem& problem) -> SmartEquilibriumResult
{
    return pimpl->solve(state, problem);
}

auto SmartEquilibriumSolver::properties() const -> const ChemicalProperties&
{
    return pimpl->properties;
}

auto SmartEquilibriumSolver::result() const -> const SmartEquilibriumResult&
{
    return pimpl->result;
}

} // namespace Reaktoro

