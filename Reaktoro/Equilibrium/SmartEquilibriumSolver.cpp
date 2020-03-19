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
#include <functional>
#include <numeric>
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
namespace detail {

/// Return the hash number of a vector.
/// https://stackoverflow.com/questions/20511347/a-good-hash-function-for-a-vector
template<typename Vec>
auto hash(Vec& vec) -> std::size_t
{
    using T = decltype(vec[0]);
    const std::hash<T> hasher;
    std::size_t seed = vec.size();
    for(const auto& i : vec)
        seed ^= hasher(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
}

} // namespace detail


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

    /// The chemical potentials at the calculated equilibrium state
    ChemicalVector u;

    /// The amounts of the species in the chemical system
    Vector n;

    /// The amounts of the equilibrium species
    Vector ne;

    /// The amounts of the elements in the equilibrium partition
    Vector be;

    /// The storage for matrix du/db = du/dn * dn/db
    Matrix dudb;

    /// The storage for vector u(iprimary)
    Vector up;

    /// The storage for matrix Mbe = inv(u(iprimary)) * du(iprimary)/db.
    Matrix Mbe;

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
        Matrix Mbe;
    };

    /// The cluster storing learned input-output data with same classification.
    struct Cluster
    {
        /// The indices of the primary species for this cluster.
        VectorXi iprimary;

        /// The hash of the indices of the primary species for this cluster.
        std::size_t label;

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

    /// Construct an SmartEquilibriumSolver::Impl instance with given chemical system.
    Impl(const ChemicalSystem& system)
    : Impl(Partition(system))
    {
    }

    /// Construct an SmartEquilibriumSolver::Impl instance with given partition of the chemical system.
    Impl(const Partition& partition)
    : system(partition.system()), partition(partition),
      properties(partition.system()), solver(partition)
    {
        // Initialize the canonicalizer with the formula matrix Ae of the equilibrium species
        canonicalizer.compute(partition.formulaMatrixEquilibriumPartition());
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

    /// Learn how to perform a full equilibrium calculation (with tracking)
    auto learn(ChemicalState& state, double T, double P, VectorConstRef be) -> void
    {
        //---------------------------------------------------------------------
        // GIBBS ENERGY MINIMIZATION CALCULATION DURING THE LEARNING PROCESS
        //---------------------------------------------------------------------
        tic(EQUILIBRIUM_STEP);

        // Perform a full chemical equilibrium calculation
        result.learning.gibbs_energy_minimization = solver.solve(state, T, P, be);

        result.timing.learning_gibbs_energy_minimization = toc(EQUILIBRIUM_STEP);

        //---------------------------------------------------------------------
        // CHEMICAL PROPERTIES UPDATE STEP DURING THE LEARNING PROCESS
        //---------------------------------------------------------------------
        tic(CHEMICAL_PROPERTIES_STEP);

        properties = solver.properties();

        result.timing.learning_chemical_properties = toc(CHEMICAL_PROPERTIES_STEP);

        //---------------------------------------------------------------------
        // SENSITIVITY MATRIX COMPUTATION STEP DURING THE LEARNING PROCESS
        //---------------------------------------------------------------------
        tic(SENSITIVITY_STEP);

        const auto& sensitivity = solver.sensitivity();

        result.timing.learning_sensitivity_matrix = toc(SENSITIVITY_STEP);

        //---------------------------------------------------------------------
        // ERROR CONTROL MATRICES ASSEMBLING STEP DURING THE LEARNING PROCESS
        //---------------------------------------------------------------------
        tic(ERROR_CONTROL_MATRICES);

        // The indices of the equilibrium species and elements
        const auto& ies = partition.indicesEquilibriumSpecies();
        const auto& iee = partition.indicesEquilibriumElements();

        // The number of equilibrium species and elements
        const auto& Ne = partition.numEquilibriumSpecies();
        const auto& Ee = partition.numEquilibriumElements();

        // The amounts of the species at the calculated equilibrium state
        n = state.speciesAmounts();

        // The amounts of the equilibrium species amounts at the calculated equilibrium state
        ne = n(ies);

        // Update the canonical form of formula matrix Ae so that we can identify primary species
        canonicalizer.updateWithPriorityWeights(ne);

        // The order of the equilibrium species as (primary, secondary)
        const auto& iorder = canonicalizer.Q();

        // Assemble the vector of indices of equilibrium species as (primary, secondary)
        VectorXi ips(ies.size());
        for(auto i = 0; i < ips.size(); ++i)  // TODO: type Indices should be alias to VectorXi, to avoid such kind of codes
            ips[i] = ies[iorder[i]];

        // The number of primary species among the equilibrium species
        const auto& Np = canonicalizer.numBasicVariables();

        // Store the indices of primary and secondary species in state
        state.equilibrium().setIndicesEquilibriumSpecies(ips, Np);

        // The indices of the primary species at the calculated equilibrium state
        const auto& iprimary = state.equilibrium().indicesPrimarySpecies();

        // The chemical potentials at the calculated equilibrium state
        u = properties.chemicalPotentials();

        // Auxiliary references to the derivatives dn/db and du/dn
        const auto& dndb = solver.sensitivity().dndb;
        const auto& dudn = u.ddn;

        // Compute the matrix du/db = du/dn * dn/db
        dudb = dudn * dndb;

        // The vector u(iprimary) with chemical potentials of primary species
        up.noalias() = u.val(iprimary);

        // The matrix du(iprimary)/dbe with derivatives of chemical potentials (primary species only)
        const auto dupdbe = dudb(iprimary, iee);

        // Compute matrix Mbe = 1/up * dup/db
        Mbe.noalias() = diag(inv(up)) * dupdbe;

        result.timing.learning_error_control_matrices = toc(ERROR_CONTROL_MATRICES);

        //---------------------------------------------------------------------
        // STORAGE STEP DURING THE LEARNING PROCESS
        //---------------------------------------------------------------------
        tic(STORAGE_STEP);

        // Generate the hash number for the indices of primary species in the state
        const auto label = detail::hash(iprimary);

        // Find the index of the cluster that has same primary species
        auto iter = std::find_if(database.clusters.begin(), database.clusters.end(),
            [&](const Cluster& cluster) { return cluster.label == label; });

        // If cluster found, store the new record in it, otherwise, create a new cluster
        if(iter < database.clusters.end())
        {
            auto& cluster = *iter;
            cluster.records.push_back({T, P, be, state, properties, sensitivity, Mbe});
            cluster.priority.extend();
        }
        else
        {
            // Create a new cluster
            Cluster cluster;
            cluster.iprimary = iprimary;
            cluster.label = label;
            cluster.records.push_back({T, P, be, state, properties, sensitivity, Mbe});
            cluster.priority.extend();

            // Append the new cluster in the database
            database.clusters.push_back(cluster);
            database.connectivity.extend();
            database.priority.extend();
        }

        result.timing.learning_storage = toc(STORAGE_STEP);
    }

    /// Estimate the equilibrium state using sensitivity derivatives (profiling the expences)
    auto estimate(ChemicalState& state, double T, double P, VectorConstRef be) -> void
    {
        // Set the estimate status to false at the beginning
        result.estimate.accepted = false;

        // Skip estimation if no cluster exists yet
        if(database.clusters.empty())
            return;

        // Auxiliary relative and absolute tolerance parameters
        const auto reltol = options.reltol;

        // The threshold used to determine elements with insignificant amounts
        const auto eps_b = options.amount_fraction_cutoff * sum(be);

        // The current set of primary species in the chemical state
        const auto& iprimary = state.equilibrium().indicesPrimarySpecies();

        Vector dbe;

        // The function that checks if a record in the database pass the error test.
        // It returns (`success`, `error`, `iprimaryspecies`), where
        // * `success` is true if error test succeeds, false otherwise.
        // * `error` is the first error value violating the tolerance
        // * `iprimaryspecies` is the index of the primary species that fails the error test
        auto pass_error_test = [&be, &dbe, &reltol, &eps_b](const Record& record) -> std::tuple<bool, double, Index>
        {
            using std::abs;
            using std::max;
            const auto& state0 = record.state;
            const auto& be0 = record.be;
            const auto& Mbe0 = record.Mbe;
            const auto& isue0 = state0.equilibrium().indicesStrictlyUnstableElements();

            dbe.noalias() = be - be0;

            // Check if state0 has strictly unstable elements (i.e. elements with zero amounts)
            // which cannot be used for Taylor estimation if positive amounts for those elements are given.
            if((dbe(isue0).array() > eps_b).any())
            {
                assert((be0(isue0).array() < eps_b).any()); // ensure this condition is never broken (effective during debug only)
                return { false, 9999, -1 };
            }

            double error = 0.0;
            const auto size = Mbe0.rows();
            for(auto i = 1; i <= size; ++i) {
                error = max(error, abs(Mbe0.row(size - i) * dbe)); // start checking primary species with least amount first
                if(error >= reltol)
                    return { false, error, size - i };
            }

            return { true, error, -1 };
        };

        // Generate the hash number for the indices of primary species in the state
        const auto label = detail::hash(iprimary);

        // The function that identifies the starting cluster index
        auto index_starting_cluster = [&]() -> Index
        {
            // If no primary species, then return number of clusters to trigger use of total usage counts of clusters
            if(iprimary.size() == 0)
                return database.clusters.size();

            // Find the index of the cluster with same set of primary species (search those with highest count first)
            for(auto icluster : database.priority.order())
                if(database.clusters[icluster].label == label)
                    return icluster;

            return database.clusters.size();
        };

        // The index of the starting cluster
        const auto icluster = index_starting_cluster();

        // The ordering of the clusters to look for (starting with icluster)
        const auto& clusters_ordering = database.connectivity.order(icluster);

        //---------------------------------------------------------------------
        // SEARCH STEP DURING THE ESTIMATE PROCESS
        //---------------------------------------------------------------------
        tic(SEARCH_STEP);

        // Iterate over all clusters starting with icluster
        for(auto jcluster : clusters_ordering)
        {
            const auto& records = database.clusters[jcluster].records;
            const auto& records_ordering = database.clusters[jcluster].priority.order();

            // Iterate over all records in current cluster
            for(auto irecord : records_ordering)
            {
                const auto& record = records[irecord];

                //---------------------------------------------------------------------
                // ERROR CONTROL STEP DURING THE ESTIMATE PROCESS
                //---------------------------------------------------------------------
                tic(ERROR_CONTROL_STEP);

                const auto [success, error, iprimaryspecies] = pass_error_test(record);

                result.timing.estimate_error_control += toc(ERROR_CONTROL_STEP);

                if(success)
                {
                    //---------------------------------------------------------------------
                    // TAYLOR PREDICTION STEP DURING THE ESTIMATE PROCESS
                    //---------------------------------------------------------------------
                    tic(TAYLOR_STEP);

                    const auto& ies = partition.indicesEquilibriumSpecies();
                    const auto& iee = partition.indicesEquilibriumElements();
                    const auto& be0 = record.be;
                    const auto& n0 = record.state.speciesAmounts();
                    const auto& y0 = record.state.equilibrium().elementChemicalPotentials();
                    const auto& z0 = record.state.equilibrium().speciesStabilities();
                    const auto& ips0 = record.state.equilibrium().indicesEquilibriumSpecies();
                    const auto& Np0 = record.state.equilibrium().numPrimarySpecies();
                    const auto& dndb0 = record.sensitivity.dndb;
                    const auto& dnedbe0 = dndb0(ies, iee);
                    const auto& ne0 = n0(ies);

                    ne.noalias() = ne0 + dnedbe0 * (be - be0);

                    n(ies) = ne;

                    const double ne_min = min(ne);
                    const double ne_sum = sum(ne);

                    const auto eps_n = options.amount_fraction_cutoff * ne_sum;

                    // Check if all projected species amounts are positive
                    if(ne_min <= -eps_n)
                        continue;

                    result.timing.estimate_search = toc(SEARCH_STEP);

                    // Set the chemical state result with estimated amounts
                    state = record.state; // ATTENTION: If this changes one day, make sure indices of equilibrium primary/secondary species, and indices of strictly unstable species/elements are also transfered from reference state to new state
                    state.setSpeciesAmounts(n);

                    // Update the chemical properties of the system
                    properties = record.properties;  // TODO: We need to estimate properties = properties0 + variation : THIS IS A TEMPORARY SOLUTION!!!

                    result.timing.estimate_taylor = toc(TAYLOR_STEP);

                    //---------------------------------------------------------------------
                    // DATABASE PRIORITY UPDATE STEP DURING THE ESTIMATE PROCESS
                    //---------------------------------------------------------------------
                    tic(PRIORITY_UPDATE_STEP);

                    // Increment priority of the current record (irecord) in the current cluster (jcluster)
                    database.clusters[jcluster].priority.increment(irecord);

                    // Increment priority of the current cluster (jcluster) with respect to starting cluster (icluster)
                    database.connectivity.increment(icluster, jcluster);

                    // Increment priority of the current cluster (jcluster)
                    database.priority.increment(jcluster);

                    result.timing.estimate_database_priority_update = toc(PRIORITY_UPDATE_STEP);

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
        const auto T = problem.temperature();
        const auto P = problem.pressure();
        const auto& iee = partition.indicesEquilibriumElements();
        be = problem.elementAmounts()(iee);
        return solve(state, T, P, be);
    }

    auto solve(ChemicalState& state, double T, double P, VectorConstRef be) -> SmartEquilibriumResult
    {
        tic(SOLVE_STEP);

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

        result.timing.solve = toc(SOLVE_STEP);

        return result;
    }

    auto solve(ChemicalState& state, double T, double P, VectorConstRef be, Index istep, Index icell) -> SmartEquilibriumResult
    {
        tic(SOLVE_STEP);

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

        result.timing.solve = toc(SOLVE_STEP);

        // Print
        if(istep == 0 || istep == 19 || istep == 2399)
        {
            const auto Ae = partition.formulaMatrixEquilibriumPartition();
            const auto& ies = partition.indicesEquilibriumSpecies();
            const auto& n = state.speciesAmounts();
            const auto& ne = n(ies);

            Vector res = abs(Ae*ne - be);
            if(icell == 0) std::cout << "res on istep " << istep << "\n" << tr(res) << std::endl;
            else std::cout << tr(res)  << std::endl;
        }
        return result;
    }

    auto outputClusterInfo() const -> void
    {
        //std::cout << "***********************************************************************************" << std::endl;
        //std::cout << "Clusters ordered by order of their creation" << std::endl;
        //std::cout << "***********************************************************************************" << std::endl;

        Index i = 0;
        for(auto cluster : database.clusters)
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
            std::cout << std::endl;
            i++;
        }
        /*
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
        */
    }
};

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

auto SmartEquilibriumSolver::solve(ChemicalState& state, double T, double P, VectorConstRef be) -> SmartEquilibriumResult
{
    return pimpl->solve(state, T, P, be);
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, double T, double P, VectorConstRef be, Index istep, Index icell) -> SmartEquilibriumResult
{
    return pimpl->solve(state, T, P, be, istep, icell);
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

