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

#include "SmartEquilibriumSolverClustering.hpp"

// C++ includes
#include <algorithm>
#include <numeric>
#include <tuple>

namespace Reaktoro {

namespace detail {

/// Return the hash number of a vector.
/// https://stackoverflow.com/questions/20511347/a-good-hash-function-for-a-vector
template<typename Vec>
auto hash(Vec &vec) -> std::size_t 
{
    using T = decltype(vec[0]);
    std::hash<T> hasher;
    std::size_t seed = vec.size();
    for (const auto &i : vec)
        seed ^= hasher(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;

}
}

SmartEquilibriumSolverClustering::SmartEquilibriumSolverClustering(const ChemicalSystem &system)
: SmartEquilibriumSolverBase(Partition(system))
{
}

SmartEquilibriumSolverClustering::SmartEquilibriumSolverClustering(const Partition& partition)
: SmartEquilibriumSolverBase(partition)
{
}

SmartEquilibriumSolverClustering::~SmartEquilibriumSolverClustering()
{
}

auto SmartEquilibriumSolverClustering::learn(ChemicalState& state, double T, double P, VectorConstRef be) -> void
{
    //---------------------------------------------------------------------
    // GIBBS ENERGY MINIMIZATION CALCULATION DURING THE LEARNING PROCESS
    //---------------------------------------------------------------------
    tic(EQUILIBRIUM_STEP)

    // Perform a full chemical equilibrium calculation
    _result.learning.gibbs_energy_minimization = solver.solve(state, T, P, be);

    // Check if the EquilibriumSolver calculation failed, if so, use cold-start
    if(!_result.learning.gibbs_energy_minimization.optimum.succeeded)
    {
        state.setSpeciesAmounts(0.0);
        _result.learning.gibbs_energy_minimization = solver.solve(state, T, P, be);
        if(!_result.learning.gibbs_energy_minimization.optimum.succeeded)
                return;
        
    }
    _result.timing.learn_gibbs_energy_minimization = toc(EQUILIBRIUM_STEP);

    //---------------------------------------------------------------------
    // CHEMICAL PROPERTIES UPDATE STEP DURING THE LEARNING PROCESS
    //---------------------------------------------------------------------
    tic(CHEMICAL_PROPERTIES_STEP)

    _properties =  solver.properties();

    _result.timing.learn_chemical_properties = toc(CHEMICAL_PROPERTIES_STEP);

    //---------------------------------------------------------------------
    // SENSITIVITY MATRIX COMPUTATION STEP DURING THE LEARNING PROCESS
    //---------------------------------------------------------------------
    tic(SENSITIVITY_STEP)

    const auto& sensitivity = solver.sensitivity();

    _result.timing.learn_sensitivity_matrix = toc(SENSITIVITY_STEP);

    //---------------------------------------------------------------------
    // ERROR CONTROL MATRICES ASSEMBLING STEP DURING THE LEARNING PROCESS
    //---------------------------------------------------------------------
    tic(ERROR_CONTROL_MATRICES)

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

    // The number of primary species among the equilibrium species (Np <= Ne)
    const auto& Np = canonicalizer.numBasicVariables();

    // Store the indices of primary and secondary species in state
    state.equilibrium().setIndicesEquilibriumSpecies(ips, Np);

    // The indices of the primary species at the calculated equilibrium state
    VectorXiConstRef iprimary = ips.head(Np);

    // The chemical potentials at the calculated equilibrium state
    u = _properties.chemicalPotentials();

    // Auxiliary references to the derivatives dn/db, dn/dT, dn/dP, and du/dn
    const auto& dndT = sensitivity.dndT;
    const auto& dndP = sensitivity.dndP;
    const auto& dndb = sensitivity.dndb;

    // Compute the matrices du/dT, du/dP, du/db
    dudT = u.ddn * dndT + u.ddT; // du/dT = ∂u/∂n*∂n/∂T + ∂u/∂T
    dudP = u.ddn * dndP + u.ddP; // du/dP = ∂u/∂n*∂n/∂P + ∂u/∂P
    dudb = u.ddn * dndb;         // du/du = ∂u/∂n*∂n/∂b

    // The vector u(iprimary) with chemical potentials of primary species
    up.noalias() = u.val(iprimary);

    // The matrix du(iprimary)/dbe with derivatives of chemical potentials (primary species only)
    const auto dupdbe = dudb(iprimary, iee);
    const auto dupdT = dudT(iprimary, 0);
    const auto dupdP = dudP(iprimary, 0);
    
    // Compute matrix Mbe = 1/up * dup/db
    Mbe.noalias() = diag(inv(up)) * dupdbe;
    MT.noalias()  = diag(inv(up)) * dupdT;
    MP.noalias()  = diag(inv(up)) * dupdP;

    _result.timing.learn_error_control_matrices = toc(ERROR_CONTROL_MATRICES);

    //---------------------------------------------------------------------
    // STORAGE STEP DURING THE LEARNING PROCESS
    //---------------------------------------------------------------------
    tic(STORAGE_STEP)

    // Generate the hash number for the indices of primary species in the state
    const auto label = detail::hash(iprimary);

    // Find the index of the cluster that has same primary species
    auto iter = std::find_if(database.clusters.begin(), database.clusters.end(),
        [&](const Cluster& cluster) { return cluster.label == label; });

    // If cluster is found, store the new record in it, otherwise, create a new cluster
    if(iter < database.clusters.end())
    {
        auto& cluster = *iter;
        cluster.records.push_back({T, P, be, state, _properties, sensitivity, Mbe, MT, MP});
        cluster.priority.extend();
    }
    else
    {
        // Create a new cluster
        Cluster cluster;
        cluster.iprimary = iprimary;
        cluster.label = label;
        cluster.records.push_back({T, P, be, state, _properties, sensitivity, Mbe, MT, MP});
        cluster.priority.extend();

        // Append the new cluster in the database
        database.clusters.push_back(cluster);
        database.connectivity.extend();
        database.priority.extend();
    }

    _result.timing.learn_storage = toc(STORAGE_STEP);
}

auto SmartEquilibriumSolverClustering::estimate(ChemicalState& state, double T, double P, VectorConstRef be) -> void
{
    // Set the estimate status to false at the beginning
    _result.estimate.accepted = false;

    // Skip estimation if no cluster exists yet
    if(database.clusters.empty())
        return;

    // Auxiliary relative and absolute tolerance parameters
    const auto reltol = options.reltol;

    // The threshold used to determine elements with insignificant amounts
    const auto eps_b = options.amount_fraction_cutoff * sum(be);

    // The current set of primary species in the chemical state
    const auto& iprimary = state.equilibrium().indicesPrimarySpecies();

    // Variations in be, T, and P
    Vector dbe;
    double dT, dP;

    // The function that checks if a record in the database pass the error test.
    // It returns (`success`, `error`, `iprimaryspecies`), where
    // * `success` is true if error test succeeds, false otherwise.
    // * `error` is the first error value violating the tolerance
    // * `iprimaryspecies` is the index of the primary species that fails the error test
    auto pass_error_test = [&be, &dbe, &T, &dT, &P, &dP, &reltol, &eps_b](const Record& record) -> std::tuple<bool, double, Index>
    {
        using std::abs;
        using std::max;
        const auto& state0 = record.state;
        const auto& be0 = record.be;
        const auto& T0 = record.T;
        const auto& P0 = record.P;
        const auto& Mbe0 = record.Mbe;
        const auto& MT0 = record.MT;
        const auto& MP0 = record.MP;
        const auto& isue0 = state0.equilibrium().indicesStrictlyUnstableElements();

        dbe.noalias() = be - be0;
        dT = T - T0;
        dP = P - P0;

        // Check if state0 has strictly unstable elements (i.e. elements with zero amounts)
        // which cannot be used for Taylor estimation if positive amounts for those elements are given.
        if((dbe(isue0).array() > eps_b).any())
        {
            assert((be0(isue0).array() < eps_b).any()); // ensure this condition is never broken (effective during debug only)
            return { false, 9999, -1 };
        }

        double error = 0.0;
        const auto size = Mbe0.rows();
        for(auto i = 1; i <= size; ++i)
        {
            const double delta_mu = Mbe0.row(size - i).dot(dbe) + MT0[size - i] * dT + MP0[size - i] * dP;
            error = max(error, abs(delta_mu)); // start checking primary species with least amount first
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

        // Find the index of the cluster with the same set of primary species (search those with highest count first)
        for(auto icluster : database.priority.order())
            if(database.clusters[icluster].label == label)
                return icluster;

        // In no cluster with the same set of primary species if found, then return number of clusters
        return database.clusters.size();
    };

    // The index of the starting cluster
    const auto icluster = index_starting_cluster();

    // The ordering of the clusters to look for (starting with icluster)
    const auto& clusters_ordering = database.connectivity.order(icluster);

    //---------------------------------------------------------------------
    // SEARCH STEP DURING THE ESTIMATE PROCESS
    //---------------------------------------------------------------------
    tic(SEARCH_STEP)

    // Iterate over all clusters (starting with icluster)
    for(auto jcluster : clusters_ordering)
    {
        // Fetch records from the cluster and the order they have to be processed in
        const auto& records = database.clusters[jcluster].records;
        const auto& records_ordering = database.clusters[jcluster].priority.order();

        // Iterate over all records in current cluster (using the  order based on the priorities)
        for(auto irecord : records_ordering)
        {
            const auto& record = records[irecord];

            //---------------------------------------------------------------------
            // ERROR CONTROL STEP DURING THE ESTIMATE PROCESS
            //---------------------------------------------------------------------
            tic(ERROR_CONTROL_STEP)

            // Check if the current record passes the error test
            const auto [success, error, iprimaryspecies] = pass_error_test(record);

            _result.timing.estimate_error_control += toc(ERROR_CONTROL_STEP);

            if(success)
            {
                _result.timing.estimate_error_control = toc(ERROR_CONTROL_STEP);

                //---------------------------------------------------------------------
                // TAYLOR PREDICTION STEP DURING THE ESTIMATE PROCESS
                //---------------------------------------------------------------------
                tic(TAYLOR_STEP)

                // Fetch reference values
                const auto& be0 = record.be;
                const auto& n0 = record.state.speciesAmounts();
                const auto& P0 = record.state.pressure();
                const auto& T0 = record.state.temperature();
                const auto& y0 = record.state.equilibrium().elementChemicalPotentials();
                const auto& z0 = record.state.equilibrium().speciesStabilities();
                const auto& ips0 = record.state.equilibrium().indicesEquilibriumSpecies();
                const auto& Np0 = record.state.equilibrium().numPrimarySpecies();
                const auto& dndb0 = record.sensitivity.dndb;
                const auto& dndP0 = record.sensitivity.dndP;
                const auto& dndT0 = record.sensitivity.dndT;

                // Fetch reference values restricted to equilibrium species only
                const auto& dnedbe0 = dndb0(ies, iee);
                const auto& dnedP0 = dndP0(ies);
                const auto& dnedT0 = dndT0(ies);
                const auto& ne0 = n0(ies);

                // Perform Taylor extrapolation
                //ne.noalias() = ne0 + dnedbe0 * (be - be0);
                ne.noalias() = ne0 + dnedbe0*(be - be0) + dnedbPe0*(P - P0) + dnedbTe0*(T - T0);

                _result.timing.estimate_taylor = toc(TAYLOR_STEP);

                // Check if all projected species amounts are positive
                const double ne_min = min(ne);
                const double ne_sum = sum(ne);
                const auto eps_n = options.amount_fraction_cutoff * ne_sum;
                if(ne_min <= -eps_n)
                    continue;

                _result.timing.estimate_search = toc(SEARCH_STEP);

                // After the search is finished successfully
                //---------------------------------------------------------------------

                // Assign small values to all the amount in the interval [cutoff, 0] (instead of mirroring above)
                for(unsigned int i = 0; i < ne.size(); ++i) 
                    if(ne[i] < 0) 
                        ne[i] = options.learning.epsilon;

                // Update the amounts of elements for the equilibrium species
                n(ies) = ne;

                // Update the chemical state res with estimated amounts
                //state = record.state; // ATTENTION: If this changes one day, make sure indices of equilibrium primary/secondary species, and indices of strictly unstable species/elements are also transfered from reference state to new state
                state.setSpeciesAmounts(ne, ies);

                // Make sure that pressure and temperature is set to the current one we are trying to predict
                state.setPressure(P);
                state.setTemperature(T);

                // Make sure that pressure and temperature is set to the current one
                state.setTemperature(T);
                state.setPressure(P);
                    
                // Update the chemical properties of the system
                _properties =  record.properties;  // TODO: We need to estimate properties = properties0 + variation : THIS IS A TEMPORARY SOLUTION!!!

                //---------------------------------------------------------------------
                // DATABASE PRIORITY UPDATE STEP DURING THE ESTIMATE PROCESS
                //---------------------------------------------------------------------
                tic(PRIORITY_UPDATE_STEP)

                // Increment priority of the current record (irecord) in the current cluster (jcluster)
                database.clusters[jcluster].priority.increment(irecord);

                // Increment priority of the current cluster (jcluster) with respect to starting cluster (icluster)
                database.connectivity.increment(icluster, jcluster);

                // Increment priority of the current cluster (jcluster)
                database.priority.increment(jcluster);

                _result.timing.estimate_database_priority_update = toc(PRIORITY_UPDATE_STEP);

                // Mark the estimated state as accepted
                _result.estimate.accepted = true;

                return;
            }
        }
    }

    _result.estimate.accepted = false;

}

auto SmartEquilibriumSolverClustering::outputInfo() const -> void
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
        //std::cout << "  RANK OF RECORDS: ";
        //for(auto j : database.clusters[i].priority.order())
        //    std::cout << database.clusters[i].priority.priorities()[j] << " ";
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

} // namespace Reaktoro

