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

#include "SmartKineticSolverClustering.hpp"

// C++ includes
#include <functional>

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
        seed ^= hasher(i) + 0x9e3779b9 + (seed << 6u) + (seed >> 2u);
    return seed;
}

/// Sort partially sorted array
void insertion_sort_from_assigned_start_(VectorXi &array, Vector values, int start)
{
    // We obtain:
    // `array` is array of indices (of species) that must be correctly shuffled in correspondence to the `values`
    // `values` is array of values (amount of species) corresponding to `indices` (in chemical state)
    // `start` is the index, from which `values` is not sorted and `array` must be shuffled in correspondence
    //
    // Example of algorithm functionality:
    //
    // array: 35      16     9    28   1    4   18  24   31          34      33
    // values 35701.7 4947.9 59.8 57.9 14.3 1.2 0.8 0.03 8.91465e-05 1459.75 560.257
    // start: 9
    //
    // 1) array[9] = 34 must be inserted between 16 and 9 because 4947.9 > 1459.75 > 59.8
    //
    // array: 35      16     34       9    28   1    4   18  24   31          33
    // values 35701.7 4947.9 1459.75  59.8 57.9 14.3 1.2 0.8 0.03 8.91465e-05 560.257
    // start : 10
    //
    // 2) array[10] = 33 must be inserted between 34 and 9 because 1459.75 > 560.257 > 59.8

    Index i_value; // the key value that must be inserted
    double amount_i_value;

    // Run starting from the first "element to be inserted" indexed by `start` till the last
    for(int i = start, head = 0; i < array.size() && head < start; head++)
    {
        // Save `array[i]` that must be inserted later in the correct place
        i_value = array[i];
        amount_i_value  = values[i];

        // If we have reached the point, where elem value can be inserted
        // i.e., the value corresponding to `i`, `values[i]`, is bigger than `values[head]`
        if(values[i] > values[head])
        {
            // Move up smaller values
            for(auto k = i; k > head; k--)
            {
                array[k] = array[k - 1];
                values[k] = values[k - 1];
            }
            // Insert elem into its proper location
            array[head] = i_value;
            values[head] = amount_i_value;

            // Change the position of the next inserted element
            i++;
            // Refresh the head index (in case next value, values[i++] is bigger)
            head = 0;
        }
        else
        {
            // Otherwise, continue marching along 'array' increasing head
            head++;
            // If `head` reached `i`, we haven't found a place to insert
            if(head == i)
            {
                head = 0;   // reset head index
                i++;        // move on to the next element to insert
            }
        }
        // Otherwise, continue marching along 'array' increasing head
    }
    return;
}

} // namespace detail

/// Implementation of the SmartKineticSolverClustering functionality

SmartKineticSolverClustering::SmartKineticSolverClustering(const ReactionSystem& reactions, const Partition& partition)
: SmartKineticSolverBase(reactions, partition)
{
}

SmartKineticSolverClustering::~SmartKineticSolverClustering()
{
}

auto SmartKineticSolverClustering::learn(ChemicalState& state, double& t, double dt) -> void
{
    //---------------------------------------------------------------------
    // INITIALIZATION STEP
    //---------------------------------------------------------------------

    // Initialize sensitivity matrix by the identity matrix
    benk_S.setIdentity();

    // Initialize the kinetics state with the data at times t0 and t0
    ODEState ode_state;
    ode_state.u0 = benk;
    ode_state.t0 = t;

    //---------------------------------------------------------------------
    // CONVENTIONAL TIME-INTEGRATION DURING THE LEARNING PROCESS
    //---------------------------------------------------------------------
    tic(INTEGRATE_STEP)

    // Calculate `benk` by the conventional numerical integration
    ode.solve(t, dt, benk, benk_S);

    _result.timing.learn_integration = toc(INTEGRATE_STEP);

    // Save the sensitivity values, the result time, and the obtain species' amount
    ode_state.t = t;
    ode_state.u = benk;
    ode_state.dudu0 = benk_S;

    // Update the composition of the amounts of equilibrium elements and kinetic species
    be = benk.head(Ee);
    nk = benk.tail(Nk);
    state.setSpeciesAmounts(nk, iks);

    //---------------------------------------------------------------------
    // ERROR CONTROL MATRICES ASSEMBLING STEP DURING THE LEARNING PROCESS
    //---------------------------------------------------------------------
    tic(ERROR_CONTROL_MATRICES)

    // The amounts of the species at the calculated equilibrium state
    Vector n = state.speciesAmounts();

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
    ChemicalVector u = _properties.chemicalPotentials();

    // Auxiliary references to the derivatives dn/db and du/dn
    const auto& dndb = equilibrium.sensitivity().dndb;
    const auto& dudn = u.ddn;

    // Compute the matrix du/db = du/dn * dn/db
    dudb = dudn * dndb;

    // The vector u(iprimary) with chemical potentials of primary species
    up.noalias() = u.val(iprimary);

    // The matrix du(iprimary)/dbe with derivatives of chemical potentials (primary species only)
    const auto dupdbe = dudb(iprimary, iee);

    // Compute matrix Mbe = 1/up * dup/db
    Mbe.noalias() = diag(inv(up)) * dupdbe;

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
        cluster.records.push_back({T, P, be, state, _properties, sensitivity, Mbe, ode_state, rates});
        cluster.priority.extend();
    }
    else
    {
        // Create a new cluster
        Cluster cluster;
        cluster.iprimary = iprimary;
        cluster.label = label;
        cluster.records.push_back({T, P, be, state, _properties, sensitivity, Mbe, ode_state, rates});
        cluster.priority.extend();

        // Append the new cluster in the database
        database.clusters.push_back(cluster);
        database.connectivity.extend();
        database.priority.extend();
    }

    _result.timing.learn_storage = toc(STORAGE_STEP);
}

auto SmartKineticSolverClustering::estimate(ChemicalState& state, double& t, double dt) -> void
{
    _result.estimate.accepted = false;

    // Skip estimation if no previous full computation has been done
    if(database.clusters.empty())
        return;

    // Relative and absolute tolerance parameters
    const auto reltol = options.tol;

    // Define initial state of the problem
    Vector benk0 = benk;
    // Define initial state of the problem
    Vector be = benk.head(Ee);
    Vector nk = benk.tail(Nk);

    const double T = state.temperature();
    const double P = state.pressure();

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
    auto pass_error_test = [&be, &dbe, &reltol, &eps_b](const KineticRecord& record) -> std::tuple<bool, double, Index>
    {
        using std::abs;
        using std::max;
        const auto& state_ref = record.state;
        const auto& be_ref = record.be;
        const auto& Mbe_ref = record.Mbe;
        const auto& isue_ref = state_ref.equilibrium().indicesStrictlyUnstableElements();

        dbe.noalias() = be - be_ref;

        // Check if state_ref has strictly unstable elements (i.e. elements with zero amounts)
        // which cannot be used for Taylor estimation if positive amounts for those elements are given.
        if((dbe(isue_ref).array() > eps_b).any())
        {
            assert((be_ref(isue_ref).array() < eps_b).any()); // ensure this condition is never broken (effective during debug only)
            return { false, 9999, -1 };
        }

        double error = 0.0;
        const auto size = Mbe_ref.rows();
        for(auto i = 1; i <= size; ++i) {
            error = max(error, abs(Mbe_ref.row(size - i) * dbe)); // start checking primary species with least amount first
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
            return database.clusters.size(); // or database.clusters.size() - 1 ?

        // Find the index of the cluster with the same set of primary species (search those with highest count first)
        for(auto icluster : database.priority.order())
            if(database.clusters[icluster].label == label)
                return icluster;

        // In no cluster with the same set of primary species if found, then return number of clusters
        return database.clusters.size(); // or database.clusters.size() - 1 ?
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
            tic(EQUILIBIRUM_SPEICES_ERROR_CONTROL_STEP)

            // Check if the current record passes the error test
            const auto [success, error, iprimaryspecies] = pass_error_test(record);

            _result.timing.estimate_error_control += toc(EQUILIBIRUM_SPEICES_ERROR_CONTROL_STEP);

            if(success)
            {
                //---------------------------------------------------------------------
                // TAYLOR DURING THE ESTIMATE PROCESS
                //---------------------------------------------------------------------
                tic(TAYLOR_STEP)

                // Fetch the data stored in the reference element
                const auto& benk0_ref = record.ode_state.u0;
                const auto& benk_ref = record.ode_state.u;
                const auto& dndn0_ref = record.ode_state.dudu0;

                // Algorithm:
                // the reference state contains:
                // u0 -> benk0_ref      is the initial condition of reference vector
                // u -> benk_ref        is already calculated by integration reference vector
                // du/du0 -> dndn0_ref  is the sensitivity w.r.t. the initial condition
                // new initial values  : benk0
                // the predicted state : benk_new = benk_ref + dndn0_ref * (benk0 - benk0_ref)

                // Perform smart estimation of benk
                Vector benk_new;
                benk_new.noalias() = benk_ref + dndn0_ref * (benk0 - benk0_ref);

                _result.timing.estimate_taylor = toc(TAYLOR_STEP);

                //---------------------------------------------------------------------
                // ERROR CONTROL THE ESTIMATE PROCESS
                //---------------------------------------------------------------------
                tic(ERROR_CONTROL_STEP)

                // Fetch the be and nk unknowns from vector benk = [be; nk]
                VectorConstRef be_new = benk_new.head(Ee);
                VectorConstRef nk_new = benk_new.tail(Nk);

                // -------------------------------------------------------------------------------------------------------------
                // Check predicted negative values of equilibrium species
                // -------------------------------------------------------------------------------------------------------------

                // Define the function checking the negativity of the equilibrium species amounts
                auto pass_negative_equilibrium_species_amounts_error_test = [&](const auto& record) -> std::tuple<bool, VectorConstRef>
                {
                    // Fetch properties of the reference state
                    const auto& state_ref = record.state;
                    const auto& n_ref = state_ref.speciesAmounts();
                    const auto& P_ref = record.state.pressure();
                    const auto& T_ref = record.state.temperature();
                    const auto& dndb_ref = record.sensitivity.dndb;
                    const auto& dndP_ref = record.sensitivity.dndP;
                    const auto& dndT_ref = record.sensitivity.dndT;
                    const auto& be_ref = benk_ref.head(Ee);

                    const auto& ne_ref = n_ref(ies);
                    const auto& dnedbe_ref = dndb_ref(ies, iee);
                    const auto& dnedP_ref = dndP_ref(ies);
                    const auto& dnedT_ref = dndT_ref(ies);

                    // Define vectors of equilibrium species ne, delta(ne), delta(lnae)
                    VectorConstRef dne = dnedbe_ref * (be_new - be_ref) + dnedP_ref * (P - P_ref) + dnedT_ref * (T - T_ref);
                    ne.noalias() = ne_ref + dne;

                    // Check if all projected species amounts are positive
                    const double ne_min = min(ne);
                    const double ne_sum = sum(ne);
                    const auto eps_n = options.amount_fraction_cutoff * ne_sum;

                    return {ne_min > -eps_n, dne};
                };
                // Check if the `negative equilibrium species amounts` pass the test
                auto [is_neg_equilibrium_test_passed, dne] = pass_negative_equilibrium_species_amounts_error_test(record);
                if(!is_neg_equilibrium_test_passed)
                    continue;

                // -------------------------------------------------------------------------------------------------------------
                // Check the variations of kinetic species
                // -------------------------------------------------------------------------------------------------------------

                // Define the function checking the variations in the kinetics rate
                auto pass_kinetic_rate_variation_error_test = [&](const auto& node, VectorConstRef dne) -> bool
                {

                    const auto& rates_ref = node.rates;
                    const auto& nk_ref = benk_ref.tail(Nk);
                    const auto& properties_ref = node.properties;

                    // Initialize delta_n = [dne; dnk]
                    Vector dnk;
                    dnk.noalias() = nk_new - nk_ref;
                    Vector dn;
                    dn.resize(Nk + Ne);
                    dn(ies) << dne;
                    dn(iks) << dnk;

                    // Initialize reaction rates
                    Vector drates = rates_ref.ddn * dn;
                    //rates.val = rates_ref.val + drates;

                    // Fetch mole fractions
                    const auto& x_ref = properties_ref.moleFractions().val;
                    VectorConstRef xk_ref = x_ref(iks);

                    for(Index i = 0; i < xk_ref.size(); ++i){
                        // If the fraction is too small, skip the variational check
                        if(xk_ref[i] < options.mole_fraction_cutoff)
                            continue;
                        if(std::abs(drates.array()[i]) > options.abstol + options.reltol * std::abs(rates_ref.val.array()[i])){
                            return false;
                        }
                    }
                    return true;
                };
                // Check if the `the variations in the kinetics rate` pass the test
                const auto is_kin_rate_variation_test_passed = pass_kinetic_rate_variation_error_test(record, dne);
                if(!is_kin_rate_variation_test_passed)
                    continue;

                _result.timing.estimate_error_control += toc(ERROR_CONTROL_STEP);

                _result.timing.estimate_search = toc(SEARCH_STEP);

                // Assign small values to all the amount in the interval [cutoff, 0] (instead of mirroring above)
                for(unsigned int i = 0; i < benk_new.size(); ++i) if(benk_new[i] < 0) benk_new[i] = options.equilibrium.epsilon;

                // -------------------------------------------------------------------------------------------------------------
                // Update the solution of kinetic problem by new estimated value
                // -------------------------------------------------------------------------------------------------------------
                benk = benk_new;
                state = record.state; // ATTENTION: If this changes one day, make sure indices of equilibrium primary/secondary species, and indices of strictly unstable species/elements are also transfered from reference state to new state

                // Update the chemical properties of the system
                //properties = record.properties;  // FIXME: We actually want to estimate properties = properties0 + variation : THIS IS A TEMPORARY SOLUTION!!!

                _result.timing.estimate_taylor = toc(TAYLOR_STEP);

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

                // Mark estimated result as accepted
                _result.estimate.accepted = true;

                // Update the time
                t += dt;

                return;
            }
        }
    }

    _result.estimate.accepted = false;

}

auto SmartKineticSolverClustering::outputInfo() const -> void
{
    std::cout << "***********************************************************************************" << std::endl;
    std::cout << "Clusters ordered by order of their creation" << std::endl;
    std::cout << "***********************************************************************************" << std::endl;

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


//    else if (options.method == SmartKineticStrategy::ClusteringExtended)
//        for(auto cluster : database_extended.clusters)
//        {
//            std::cout << "CLUSTER #" << i << std::endl;
//            std::cout << "  RANK OF CLUSTER: " << database_extended.priority.priorities()[i] << std::endl;
//            std::cout << "  PRIMARY SPECIES: ";
//            for(auto j : database_extended.clusters[i].iprimary_ikin)
//                std::cout << system.species(j).name() << " ";
//            std::cout << std::endl;
//            std::cout << "  NUMBER OF RECORDS: " << database_extended.clusters[i].records.size() << std::endl;
//            //std::cout << "  RANK OF RECORDS: ";
//            //for(auto j : database.clusters[i].priority.order())
//            //    std::cout << database.clusters[i].priority.priorities()[j] << " ";
//            std::cout << std::endl;
//            std::cout << std::endl;
//            i++;
//        }
    std::cout << "***********************************************************************************" << std::endl;

    // Output content of the smart equilibrium solver if it is used
    if(options.use_smart_equilibrium_solver)
        smart_equilibrium.outputInfo();

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

/// Implementation of the SmartKineticSolverClusteringExtended functionality

SmartKineticSolverClusteringExtended::SmartKineticSolverClusteringExtended(const ReactionSystem& reactions, const Partition& partition)
: SmartKineticSolverBase(reactions, partition)
{
}

SmartKineticSolverClusteringExtended::~SmartKineticSolverClusteringExtended()
{
}

auto SmartKineticSolverClusteringExtended::learn(ChemicalState& state, double& t, double dt) -> void
{
    //---------------------------------------------------------------------
    // INITIALIZATION STEP
    //---------------------------------------------------------------------

    // Initialize sensitivity matrix by the identity matrix
    benk_S.setIdentity();

    // Initialize the kinetics state with the data at times t0 and t0
    ODEState ode_state;
    ode_state.u0 = benk;
    ode_state.t0 = t;

    //---------------------------------------------------------------------
    // CONVENTIONAL TIME-INTEGRATION DURING THE LEARNING PROCESS
    //---------------------------------------------------------------------
    tic(INTEGRATE_STEP)

    // Calculate `benk` by the conventional numerical integration
    ode.solve(t, dt, benk, benk_S);
    // Correct negative values by replacing them with small epsilon
    for(unsigned int i = 0; i < benk.size(); ++i) if(benk[i] < 0) benk[i] = options.equilibrium.epsilon;
    //ode.solve_implicit_1st_order(t, dt, benk, benk_S); // TODO: compare CVODE and 1st order implicit scheme
    //ode.integrate(t, benk, t + dt, benk_S); // TODO: produces a delayed estimations, why?
    _result.timing.learn_integration += toc(INTEGRATE_STEP);

    // Save the sensitivity values, the result time, and the obtain species' amount
    ode_state.t = t;
    ode_state.u = benk;
    ode_state.dudu0 = benk_S;

    // Update the composition of the amounts of equilibrium elements and kinetic species
    be = benk.head(Ee);
    nk = benk.tail(Nk);
    state.setSpeciesAmounts(nk, iks);

    //---------------------------------------------------------------------
    // ERROR CONTROL MATRICES ASSEMBLING STEP DURING THE LEARNING PROCESS
    //---------------------------------------------------------------------
    tic(ERROR_CONTROL_MATRICES)

    // The amounts of the species at the calculated equilibrium state
    Vector n = state.speciesAmounts();

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

    // The indices of the primary species in the calculated equilibrium state
    VectorXiConstRef iprimary = ips.head(Np);

    //---------------------------------------------------------------------
    // CONSTRUCT AN EXTENDED IPRIMARY_IKIN VECTOR FOR CLUSTERING
    //---------------------------------------------------------------------

    // The indices of the primary and kinetic species the calculated equilibrium state
    VectorXi iprimary_ikin_tmp = VectorXi::Zero(iprimary.size() + iks.size());
    for (auto i = 0; i < iprimary.size(); ++i) iprimary_ikin_tmp[i] = iprimary[i];
    for (auto i = iprimary.size(); i < iprimary_ikin_tmp.size(); ++i) iprimary_ikin_tmp[i] = iks[i - iprimary.size()];

    // The amounts of species corresponding to indices iprimary_ikin
    Vector ne_iprimary_iks = Vector::Zero(iprimary_ikin_tmp.size());
    for (auto i = 0; i < iprimary_ikin_tmp.size(); ++i)
        ne_iprimary_iks[i] = state.speciesAmount(iprimary_ikin_tmp[i]);

    // Sort iprimary_ikin_tmp in correspondence with the values ne_iprimary_iks
    detail::insertion_sort_from_assigned_start_(iprimary_ikin_tmp, ne_iprimary_iks, iprimary.size());

    // The indices of the primary and kinetic species the calculated equilibrium state (in correspondence with the values ne_iprimary_iks)
    VectorXiConstRef iprimary_ikin = iprimary_ikin_tmp;

    //---------------------------------------------------------------------

    // The chemical potentials at the calculated equilibrium state
    ChemicalVector u = _properties.chemicalPotentials();

    // Auxiliary references to the derivatives dn/db and du/dn
    const auto& dndb = equilibrium.sensitivity().dndb;
    //const auto& dndb = options.use_smart_equilibrium_solver ? smart_equilibrium.sensitivity().dndb : equilibrium.sensitivity().dndb;
    const auto& dudn = u.ddn;

    // Compute the matrix du/db = du/dn * dn/db
    dudb = dudn * dndb;

    // The vector u(iprimary) with chemical potentials of primary species
    up.noalias() = u.val(iprimary);

    // The matrix du(iprimary)/dbe with derivatives of chemical potentials (primary species only)
    const auto dupdbe = dudb(iprimary, iee);

    // Compute matrix Mbe = 1/up * dup/db
    Mbe.noalias() = diag(inv(up)) * dupdbe;

    _result.timing.learn_error_control_matrices = toc(ERROR_CONTROL_MATRICES);

    //---------------------------------------------------------------------
    // STORAGE STEP DURING THE LEARNING PROCESS
    //---------------------------------------------------------------------
    tic(STORAGE_STEP)

    // Generate the hash number for the indices of primary species in the state
    const auto label = detail::hash(iprimary);

    // Generate the hash number for the indices of primary and kinetic species in the state
    const auto label_extended = detail::hash(iprimary_ikin);

    // Find the index of the cluster that has same primary species
    auto iter = std::find_if(database_extended.clusters.begin(), database_extended.clusters.end(),
                             [&](const ClusterExtended& cluster) { return cluster.label_extended == label_extended; });

    // If cluster is found, store the new record in it, otherwise, create a new cluster
    if(iter < database_extended.clusters.end())
    {
        auto& cluster = *iter;
        //cluster.records.push_back({T, P, be0, state, properties, sensitivity, Mbe, ode_state, rates});
        cluster.records.push_back({T, P, be, state, _properties, sensitivity, Mbe, ode_state, rates});
        cluster.priority.extend();
    }
    else
    {
        // Create a new cluster
        ClusterExtended cluster;
        cluster.iprimary = iprimary;
        cluster.label = label;
        cluster.iprimary_ikin = iprimary_ikin;
        cluster.label_extended = label_extended;

        //cluster.records.push_back({T, P, be0, state, properties, sensitivity, Mbe, ode_state, rates});
        cluster.records.push_back({T, P, be, state, _properties, sensitivity, Mbe, ode_state, rates});
        cluster.priority.extend();

        // Append the new cluster in the database
        database_extended.clusters.push_back(cluster);
        database_extended.connectivity.extend();
        database_extended.priority.extend();
    }

    _result.timing.learn_storage = toc(STORAGE_STEP);
}

auto SmartKineticSolverClusteringExtended::estimate(ChemicalState& state, double& t, double dt) -> void
{
    _result.estimate.accepted = false;

    // Skip estimation if no previous full computation has been done
    if(database_extended.clusters.empty())
        return;

    // Relative and absolute tolerance parameters
    const auto reltol = options.tol;

    // Define initial state of the problem
    Vector benk0 = benk;
    // Define initial state of the problem
    Vector be = benk.head(Ee);
    Vector nk = benk.tail(Nk);

    const double T = state.temperature();
    const double P = state.pressure();

    // The threshold used to determine elements with insignificant amounts
    const auto eps_b = options.amount_fraction_cutoff * sum(be);

    // The current set of primary species in the chemical state
    const auto& iprimary = state.equilibrium().indicesPrimarySpecies();

    //---------------------------------------------------------------------
    // CONSTRUCT AN EXTENDED IPRIMARY_IKIN VECTOR FOR CLUSTERING
    //---------------------------------------------------------------------

    // The indices of the primary and kinetic species the calculated equilibrium state
    VectorXi iprimary_ikin_tmp = VectorXi::Zero(iprimary.size() + iks.size());
    for (auto i = 0; i < iprimary.size(); ++i) iprimary_ikin_tmp[i] = iprimary[i];
    for (auto i = iprimary.size(); i < iprimary_ikin_tmp.size(); ++i) iprimary_ikin_tmp[i] = iks[i - iprimary.size()];

    // The amounts of species corresponding to indices `iprimary_ikin`
    Vector ne_iprimary_iks = Vector::Zero(iprimary_ikin_tmp.size());
    for (auto i = 0; i < iprimary_ikin_tmp.size(); ++i)
        ne_iprimary_iks[i] = state.speciesAmount(iprimary_ikin_tmp[i]);

    // Sort `iprimary_ikin_tmp` in correspondence with the values `ne_iprimary_iks`
    // if `iprimary` is not empty, i.e., [35 16 9 28 18  4  1 24 31 | 34 33],
    // start sorting from the first element in the second part of the `iprimary_ikin_tmp` array
    // otherwise, [| 34 33], start from the second element
    if(iprimary.size())
        detail::insertion_sort_from_assigned_start_(iprimary_ikin_tmp, ne_iprimary_iks, iprimary.size());
    else if (iks.size())
        detail::insertion_sort_from_assigned_start_(iprimary_ikin_tmp, ne_iprimary_iks, 1);

    // The indices of the primary and kinetic species the calculated equilibrium state (in correspondence with the values ne_iprimary_iks)
    VectorXiConstRef iprimary_ikin = iprimary_ikin_tmp;
    //---------------------------------------------------------------------

    Vector dbe;

    // The function that checks if a record in the database pass the error test.
    // It returns (`success`, `error`, `iprimaryspecies`), where
    // * `success` is true if error test succeeds, false otherwise.
    // * `error` is the first error value violating the tolerance
    // * `iprimaryspecies` is the index of the primary species that fails the error test
    auto pass_error_test = [&be, &dbe, &reltol, &eps_b](const KineticRecord& record) -> std::tuple<bool, double, Index>
    {
        using std::abs;
        using std::max;
        const auto& state_ref = record.state;
        const auto& be_ref = record.be;
        const auto& Mbe_ref = record.Mbe;
        const auto& isue_ref = state_ref.equilibrium().indicesStrictlyUnstableElements();

        dbe.noalias() = be - be_ref;

        // Check if state_ref has strictly unstable elements (i.e. elements with zero amounts)
        // which cannot be used for Taylor estimation if positive amounts for those elements are given.
        if((dbe(isue_ref).array() > eps_b).any())
        {
            assert((be_ref(isue_ref).array() < eps_b).any()); // ensure this condition is never broken (effective during debug only)
            return { false, 9999, -1 };
        }

        double error = 0.0;
        const auto size = Mbe_ref.rows();
        for(auto i = 1; i <= size; ++i) {
            error = max(error, abs(Mbe_ref.row(size - i) * dbe)); // start checking primary species with least amount first
            if(error >= reltol)
                return { false, error, size - i };
        }

        return { true, error, -1 };
    };

    // Generate the hash number for the indices of primary species in the state
    const auto label_extended = detail::hash(iprimary_ikin);

    // The function that identifies the starting cluster index
    auto index_starting_cluster = [&]() -> Index
    {
        // If no primary species, then return number of clusters to trigger use of total usage counts of clusters
        if(iprimary.size() == 0)
            return database_extended.clusters.size(); // or database.clusters.size() - 1 ?

        // Find the index of the cluster with the same set of primary species (search those with highest count first)
        for(auto icluster : database_extended.priority.order())
            if(database_extended.clusters[icluster].label_extended == label_extended)
                return icluster;

        // In no cluster with the same set of primary species if found, then return number of clusters
        return database_extended.clusters.size(); // or database.clusters.size() - 1 ?
    };

    // The index of the starting cluster
    const auto icluster = index_starting_cluster();

    // The ordering of the clusters to look for (starting with icluster)
    const auto& clusters_ordering = database_extended.connectivity.order(icluster);

    //---------------------------------------------------------------------
    // SEARCH STEP DURING THE ESTIMATE PROCESS
    //---------------------------------------------------------------------
    tic(SEARCH_STEP)

    int cluster_count = 0;
    // Iterate over all clusters (starting with icluster)
    for(auto jcluster : clusters_ordering)
    {
        // Fetch records from the cluster and the order they have to be processed in
        const auto& records = database_extended.clusters[jcluster].records;
        const auto& records_ordering = database_extended.clusters[jcluster].priority.order();

        int record_count = 0;
        // Iterate over all records in current cluster (using the  order based on the priorities)
        for(auto irecord : records_ordering)
        {
            const auto& record = records[irecord];

            //---------------------------------------------------------------------
            // ERROR CONTROL STEP DURING THE ESTIMATE PROCESS
            //---------------------------------------------------------------------
            tic(EQUILIBIRUM_SPEICES_ERROR_CONTROL_STEP)

            // Check if the current record passes the error test
            const auto [success, error, iprimaryspecies] = pass_error_test(record);

            if(success)
            {
                _result.timing.estimate_error_control += toc(EQUILIBIRUM_SPEICES_ERROR_CONTROL_STEP);

                //---------------------------------------------------------------------
                // TAYLOR DURING THE ESTIMATE PROCESS
                //---------------------------------------------------------------------
                tic(TAYLOR_STEP)

                // Fetch the data stored in the reference element
                const auto& benk0_ref = record.ode_state.u0;
                const auto& benk_ref = record.ode_state.u;
                const auto& dndn0_ref = record.ode_state.dudu0;

                // Algorithm:
                // the reference state contains:
                // u0 -> benk0_ref      is the initial condition of reference vector
                // u -> benk_ref        is already calculated by integration reference vector
                // du/du0 -> dndn0_ref  is the sensitivity w.r.t. the initial condition
                // new initial values  : benk0
                // the predicted state : benk_new = benk_ref + dndn0_ref * (benk0 - benk0_ref)

                // Perform smart estimation of benk
                Vector benk_new;
                benk_new.noalias() = benk_ref + dndn0_ref * (benk0 - benk0_ref);

                _result.timing.estimate_taylor = toc(TAYLOR_STEP);

                //---------------------------------------------------------------------
                // ERROR CONTROL THE ESTIMATE PROCESS
                //---------------------------------------------------------------------
                tic(ERROR_CONTROL_STEP)

                // Fetch the be and nk unknowns from vector benk = [be; nk]
                VectorConstRef be_new = benk_new.head(Ee);
                VectorConstRef nk_new = benk_new.tail(Nk);

                // -------------------------------------------------------------------------------------------------------------
                // Check predicted negative values of equilibrium species
                // -------------------------------------------------------------------------------------------------------------

                // Define the function checking the negativity of the equilibrium species amounts
                auto pass_negative_equilibrium_species_amounts_error_test = [&](const auto& record) -> std::tuple<bool, VectorConstRef, double, double>
                {
                    // Fetch properties of the reference state
                    const auto& state_ref = record.state;
                    const auto& n_ref = state_ref.speciesAmounts();
                    const auto& P_ref = record.state.pressure();
                    const auto& T_ref = record.state.temperature();
                    const auto& dndb_ref = record.sensitivity.dndb;
                    const auto& dndP_ref = record.sensitivity.dndP;
                    const auto& dndT_ref = record.sensitivity.dndT;
                    const auto& be_ref = benk_ref.head(Ee);

                    const auto& ne_ref = n_ref(ies);
                    const auto& dnedbe_ref = dndb_ref(ies, iee);
                    const auto& dnedP_ref = dndP_ref(ies);
                    const auto& dnedT_ref = dndT_ref(ies);

                    // Define vectors of equilibrium species ne, delta(ne), delta(lnae)
                    VectorConstRef dne = dnedbe_ref * (be_new - be_ref) + dnedP_ref * (P - P_ref) + dnedT_ref * (T - T_ref);
                    ne.noalias() = ne_ref + dne;

                    // Check if all projected species amounts are positive
                    const double ne_min = min(ne);
                    const double ne_sum = sum(ne);
                    const auto eps_n = options.amount_fraction_cutoff * ne_sum;

                    return {ne_min > -eps_n, dne, ne_min, eps_n};
                    //return {ne.minCoeff() > options.cutoff, dne};
                };
                // Check if the `negative equilibrium species amounts` pass the test
                auto [is_neg_equilibrium_test_passed, dne, ne_min, eps_n] = pass_negative_equilibrium_species_amounts_error_test(record);
                if(!is_neg_equilibrium_test_passed){
                    continue;
                }

                // -------------------------------------------------------------------------------------------------------------
                // Check the variations of kinetic species
                // -------------------------------------------------------------------------------------------------------------

                // Define the function checking the variations in the kinetics rate
                auto pass_kinetic_rate_variation_error_test = [&](const auto& node, VectorConstRef dne) -> bool
                {

                    const auto& rates_ref = node.rates;
                    const auto& nk_ref = benk_ref.tail(Nk);
                    const auto& properties_ref = node.properties;

                    // Initialize delta_n = [dne; dnk]
                    Vector dnk;
                    dnk.noalias() = nk_new - nk_ref;
                    Vector dn;
                    dn.resize(Nk + Ne);
                    dn(ies) << dne;
                    dn(iks) << dnk;

                    // Initialize reaction rates
                    Vector drates = rates_ref.ddn * dn;
                    Vector drates_kin = rates_ref.ddn(Eigen::all, iks) * dnk;

                    // Fetch mole fractions
                    const auto& x_ref = properties_ref.moleFractions().val;
                    VectorConstRef xk_ref = x_ref(iks);

                    for(Index i = 0; i < xk_ref.size(); ++i){
                        // If the fraction is too small, skip the variational check
                        if(xk_ref[i] < options.mole_fraction_cutoff)
                            continue;
                        if(std::abs(drates.array()[i]) > options.abstol + options.reltol * std::abs(rates_ref.val.array()[i])){
                            return false;
                        }
                    }
                    return true;
                };
                // Check if the `the variations in the kinetics rate` pass the test
                const auto is_kin_rate_variation_test_passed = pass_kinetic_rate_variation_error_test(record, dne);
                if(!is_kin_rate_variation_test_passed){
                    continue;
                }

                _result.timing.estimate_error_control += toc(ERROR_CONTROL_STEP);

                _result.timing.estimate_search = toc(SEARCH_STEP);

                // Assign small values to all the amount in the interval [cutoff, 0] (instead of mirroring above)
                for(unsigned int i = 0; i < benk_new.size(); ++i) if(benk_new[i] < 0) benk_new[i] = options.equilibrium.epsilon;

                // -------------------------------------------------------------------------------------------------------------
                // Update the solution of kinetic problem by new estimated value
                // -------------------------------------------------------------------------------------------------------------
                benk = benk_new;
                state = record.state; // ATTENTION: If this changes one day, make sure indices of equilibrium primary/secondary species, and indices of strictly unstable species/elements are also transfered from reference state to new state

                // Update the chemical properties of the system
                //properties = record.properties;  // FIXME: We actually want to estimate properties = properties0 + variation : THIS IS A TEMPORARY SOLUTION!!!

                _result.timing.estimate_taylor = toc(TAYLOR_STEP);

                //---------------------------------------------------------------------
                // DATABASE PRIORITY UPDATE STEP DURING THE ESTIMATE PROCESS
                //---------------------------------------------------------------------
                tic(PRIORITY_UPDATE_STEP)

                // Increment priority of the current record (irecord) in the current cluster (jcluster)
                database_extended.clusters[jcluster].priority.increment(irecord);

                // Increment priority of the current cluster (jcluster) with respect to starting cluster (icluster)
                database_extended.connectivity.increment(icluster, jcluster);

                // Increment priority of the current cluster (jcluster)
                database_extended.priority.increment(jcluster);

                _result.timing.estimate_database_priority_update = toc(PRIORITY_UPDATE_STEP);

                // Mark estimated result as accepted
                _result.estimate.accepted = true;

                // Update the time
                t += dt;

                return;
            }
            record_count++;
        }
        cluster_count++;
    }

    _result.estimate.accepted = false;

}

auto SmartKineticSolverClusteringExtended::outputInfo() const -> void

{
    std::cout << "***********************************************************************************" << std::endl;
    std::cout << "Clusters ordered by order of their creation" << std::endl;
    std::cout << "***********************************************************************************" << std::endl;

    Index i = 0;
    for(auto cluster : database_extended.clusters)
    {
        std::cout << "CLUSTER #" << i << std::endl;
        std::cout << "  RANK OF CLUSTER: " << database_extended.priority.priorities()[i] << std::endl;
        std::cout << "  PRIMARY SPECIES: ";
        for(auto j : database_extended.clusters[i].iprimary)
            std::cout << system.species(j).name() << " ";
        std::cout << std::endl;
        std::cout << "  NUMBER OF RECORDS: " << database_extended.clusters[i].records.size() << std::endl;
        //std::cout << "  RANK OF RECORDS: ";
        //for(auto j : database.clusters[i].priority.order())
        //    std::cout << database.clusters[i].priority.priorities()[j] << " ";
        std::cout << std::endl;
        std::cout << std::endl;
        i++;
    }


//    else if (options.method == SmartKineticStrategy::ClusteringExtended)
//        for(auto cluster : database_extended.clusters)
//        {
//            std::cout << "CLUSTER #" << i << std::endl;
//            std::cout << "  RANK OF CLUSTER: " << database_extended.priority.priorities()[i] << std::endl;
//            std::cout << "  PRIMARY SPECIES: ";
//            for(auto j : database_extended.clusters[i].iprimary_ikin)
//                std::cout << system.species(j).name() << " ";
//            std::cout << std::endl;
//            std::cout << "  NUMBER OF RECORDS: " << database_extended.clusters[i].records.size() << std::endl;
//            //std::cout << "  RANK OF RECORDS: ";
//            //for(auto j : database.clusters[i].priority.order())
//            //    std::cout << database.clusters[i].priority.priorities()[j] << " ";
//            std::cout << std::endl;
//            std::cout << std::endl;
//            i++;
//        }
    std::cout << "***********************************************************************************" << std::endl;

    // Output content of the smart equilibrium solver if it is used
    if(options.use_smart_equilibrium_solver)
        smart_equilibrium.outputInfo();

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
