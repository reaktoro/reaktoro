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

#include "SmartEquilibriumSolver.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Profiling.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConditions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumPredictor.hpp>
#include <Reaktoro/Equilibrium/EquilibriumRestrictions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumResult.hpp>

namespace Reaktoro {
namespace detail {

/// Return the hash number of a vector.
/// https://stackoverflow.com/questions/20511347/a-good-hash-function-for-a-vector
template<typename Vec>
auto hash(Vec& vec) -> std::size_t
{
    using T = decltype(vec[0]);
    std::hash<T> hasher;
    std::size_t seed = vec.size();
    for(auto const& i : vec)
        seed ^= hasher(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
}

} // namespace detail

struct SmartEquilibriumSolver::Impl
{
    EquilibriumSolver solver;

    EquilibriumSensitivity sensitivity;

    EquilibriumConditions conditions;

    SmartEquilibriumOptions options;

    SmartEquilibriumResult result;

    /// The database with learned input-output data points.
    SmartEquilibriumSolver::Database database;

    /// Construct a SmartEquilibriumSolver::Impl object with given equilibrium problem specifications.
    Impl(EquilibriumSpecs const& specs)
    : solver(specs), sensitivity(specs), conditions(specs)
    {
        // Initialize the equilibrium solver with the default options
        setOptions(options);
    }

    //=================================================================================================================
    //
    // CHEMICAL EQUILIBRIUM METHODS
    //
    //=================================================================================================================

    auto solve(ChemicalState& state) -> SmartEquilibriumResult
    {
        conditions.temperature(state.temperature());
        conditions.pressure(state.pressure());
        conditions.surfaceAreas(state.surfaceAreas());
        return solve(state, conditions);
    }

    auto solve(ChemicalState& state, EquilibriumRestrictions const& restrictions) -> SmartEquilibriumResult
    {
        errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumRestrictions is currently not supported.");
        return {};
    }

    auto solve(ChemicalState& state, EquilibriumConditions const& conditions) -> SmartEquilibriumResult
    {
        const ArrayXd c0 = solver.componentAmounts(state);
        return solve(state, conditions, c0);
    }

    auto solve(ChemicalState& state, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions) -> SmartEquilibriumResult
    {
        errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumRestrictions is currently not supported.");
        return {};
    }

    //=================================================================================================================
    //
    // CHEMICAL EQUILIBRIUM METHODS WITH SENSITIVITY CALCULATION
    //
    //=================================================================================================================

    auto solve(ChemicalState& state, EquilibriumSensitivity& sensitivity) -> SmartEquilibriumResult
    {
        errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumSensitivity is currently not supported.");
        return {};
    }

    auto solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, EquilibriumRestrictions const& restrictions) -> SmartEquilibriumResult
    {
        errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumSensitivity is currently not supported.");
        return {};
    }

    auto solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, EquilibriumConditions const& conditions) -> SmartEquilibriumResult
    {
        errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumSensitivity is currently not supported.");
        return {};
    }

    auto solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions) -> SmartEquilibriumResult
    {
        errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumSensitivity is currently not supported.");
        return {};
    }

    //=================================================================================================================
    //
    // CHEMICAL EQUILIBRIUM METHODS WITH GIVEN AMOUNTS OF CONSERVATIVE COMPONENTS
    //
    //=================================================================================================================

    auto solve(ChemicalState& state, ArrayXdConstRef c0) -> SmartEquilibriumResult
    {
        conditions.temperature(state.temperature());
        conditions.pressure(state.pressure());
        conditions.surfaceAreas(state.surfaceAreas());
        return solve(state, conditions, c0);
    }

    auto solve(ChemicalState& state, EquilibriumRestrictions const& restrictions, ArrayXdConstRef c0) -> SmartEquilibriumResult
    {
        errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumRestrictions is currently not supported.");
        return {};
    }

    auto solve(ChemicalState& state, EquilibriumConditions const& conditions, ArrayXdConstRef c0) -> SmartEquilibriumResult
    {
        tic(SOLVE_STEP)

        // Absolutely ensure an exact Hessian of the Gibbs energy function is used in the calculations
        setOptions(options);

        // Reset the result of the last smart equilibrium calculation
        result = {};

        // Perform a smart estimate of the chemical state
        timeit( predict(state, conditions, c0), result.timing.estimate= )

        // Perform a learning step if the smart prediction is not satisfactory
        if(!result.estimate.accepted)
            timeit( learn(state, conditions, c0), result.timing.learn= )

        result.timing.solve = toc(SOLVE_STEP);

        return result;
    }

    auto solve(ChemicalState& state, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions, ArrayXdConstRef c0) -> SmartEquilibriumResult
    {
        errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumRestrictions is currently not supported.");
        return {};
    }

    //=================================================================================================================
    //
    // CHEMICAL EQUILIBRIUM METHODS WITH GIVEN AMOUNTS OF CONSERVATIVE COMPONENTS AND SENSITIVITY CALCULATION
    //
    //=================================================================================================================

    auto solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, ArrayXdConstRef c0) -> SmartEquilibriumResult
    {
        errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumSensitivity is currently not supported.");
        return {};
    }

    auto solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, EquilibriumRestrictions const& restrictions, ArrayXdConstRef c0) -> SmartEquilibriumResult
    {
        errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumSensitivity is currently not supported.");
        return {};
    }

    auto solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, EquilibriumConditions const& conditions, ArrayXdConstRef c0) -> SmartEquilibriumResult
    {
        errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumSensitivity is currently not supported.");
        return {};
    }

    auto solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions, ArrayXdConstRef c0) -> SmartEquilibriumResult
    {
        errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumSensitivity is currently not supported.");
        return {};
    }

    //=================================================================================================================
    //
    // LEARN AND PREDICT METHODS
    //
    //=================================================================================================================

    /// Perform a learning operation in which a full chemical equilibrium calculation is performed.
    auto learn(ChemicalState& state, EquilibriumConditions const& conditions, ArrayXdConstRef const& c0) -> void
    {
        //---------------------------------------------------------------------
        // GIBBS ENERGY MINIMIZATION CALCULATION DURING THE LEARNING PROCESS
        //---------------------------------------------------------------------
        tic(EQUILIBRIUM_STEP)

        // Perform a full chemical equilibrium solve with sensitivity derivatives calculation
        result.learning.gibbs_energy_minimization = solver.solve(state, sensitivity, conditions, c0);

        result.timing.learn_gibbs_energy_minimization = toc(EQUILIBRIUM_STEP);

        //---------------------------------------------------------------------
        // ERROR CONTROL MATRICES ASSEMBLING STEP DURING THE LEARNING PROCESS
        //---------------------------------------------------------------------

        // The indices of the primary species at the calculated equilibrium state
        auto const& iprimary = state.equilibrium().indicesPrimarySpecies();

        // Create an equilibrium predictor object with computed equilibrium state and its sensitivities
        EquilibriumPredictor predictor(state, sensitivity);

        //---------------------------------------------------------------------
        // STORAGE STEP DURING THE LEARNING PROCESS
        //---------------------------------------------------------------------
        tic(STORAGE_STEP)

        // Generate the hash number for the indices of primary species in the state
        const auto label = detail::hash(iprimary);

        // Find the index of the cluster that has same primary species
        // auto iter = std::find_if(database.clusters.begin(), database.clusters.end(),
        //     [&](Cluster const& cluster) { return cluster.label == label; });
        auto icluster = indexfn(database.clusters, RKT_LAMBDA(cluster, cluster.label == label));


        // If cluster is found, store the new record in it, otherwise, create a new cluster
        if(icluster < database.clusters.size())
        {
            auto& cluster = database.clusters[icluster];
            cluster.records.push_back({ state, conditions, sensitivity, predictor });
            cluster.priority.extend();
        }
        else
        {
            // Create a new cluster
            Cluster cluster;
            cluster.iprimary = iprimary;
            cluster.label = label;
            cluster.records.push_back({ state, conditions, sensitivity, predictor });
            cluster.priority.extend();

            // Append the new cluster in the database
            database.clusters.push_back(cluster);
            database.connectivity.extend();
            database.priority.extend();
        }

        result.timing.learn_storage = toc(STORAGE_STEP);
    }

    /// Perform a prediction operation in which a chemical equilibrium state is estimated using a first-order Taylor approximation.
    auto predict(ChemicalState& state, EquilibriumConditions const& conditions, ArrayXdConstRef const& c) -> void
    {
        // Set the estimate status to false at the beginning
        result.estimate.accepted = false;

        // Skip estimation if no cluster exists yet
        if(database.clusters.empty())
            return;

        // Auxiliary relative and absolute tolerance parameters
        const auto reltol = options.reltol;

        // The current set of primary species in the chemical state
        auto const& iprimary = state.equilibrium().indicesPrimarySpecies();

        // The number of primary species
        auto const numprimary = iprimary.size();

        // The function that checks if a record in the database pass the error test.
        // It returns (`success`, `error`, `iprimaryspecies`), where
        // * `success` is true if error test succeeds, false otherwise.
        // * `error` is the first error value violating the tolerance
        // * `iprimaryspecies` is the index of the primary species that fails the error test
        auto pass_error_test = [&numprimary, &iprimary, &conditions, &c, &reltol](Record const& record) -> Tuple<bool, double, Index>
        {
            using std::abs;
            using std::max;

            auto const& predictor = record.predictor;

            double error = 0.0;
            for(auto i = 1; i <= numprimary; ++i)
            {
                const auto ispecies = iprimary[numprimary - i];

                const auto mu0 = predictor.referenceSpeciesChemicalPotential(ispecies);
                const auto mu1 = predictor.predictSpeciesChemicalPotential(ispecies, conditions, c);

                const auto delta_mu = (mu0 - mu1)/mu0;

                error = max(error, abs(delta_mu)); // start checking primary species with least amount first
                if(error >= reltol)
                    return { false, error, numprimary - i };
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
        auto const& clusters_ordering = database.connectivity.order(icluster);

        //---------------------------------------------------------------------
        // SEARCH STEP DURING THE ESTIMATE PROCESS
        //---------------------------------------------------------------------
        tic(SEARCH_STEP)

        // Iterate over all clusters (starting with icluster)
        for(auto jcluster : clusters_ordering)
        {
            // Fetch records from the cluster and the order they have to be processed in
            auto const& records = database.clusters[jcluster].records;
            auto const& records_ordering = database.clusters[jcluster].priority.order();

            // Iterate over all records in current cluster (using the  order based on the priorities)
            for(auto irecord : records_ordering)
            {
                auto const& record = records[irecord];

                //---------------------------------------------------------------------
                // ERROR CONTROL STEP DURING THE ESTIMATE PROCESS
                //---------------------------------------------------------------------
                tic(ERROR_CONTROL_STEP)

                // Check if the current record passes the error test
                const auto [success, error, iprimaryspecies] = pass_error_test(record);

                result.timing.estimate_error_control += toc(ERROR_CONTROL_STEP);

                if(success)
                {
                    result.timing.estimate_error_control = toc(ERROR_CONTROL_STEP);

                    //---------------------------------------------------------------------
                    // TAYLOR PREDICTION STEP DURING THE ESTIMATE PROCESS
                    //---------------------------------------------------------------------
                    tic(TAYLOR_STEP)

                    auto const& predictor = record.predictor;

                    predictor.predict(state, conditions, c);

                    result.timing.estimate_taylor = toc(TAYLOR_STEP);

                    // Check if all projected species amounts are positive
                    auto const& n = state.speciesAmounts();
                    const double nmin = n.minCoeff();
                    const double nsum = n.sum();
                    const auto eps_n = options.amount_fraction_cutoff * nsum;
                    if(nmin <= -eps_n)
                        continue;

                    result.timing.estimate_search = toc(SEARCH_STEP);

                    // After the search is finished successfully
                    //---------------------------------------------------------------------

                    // Assign small values to all the amount in the interval [cutoff, 0] (instead of mirroring above)
                    for(auto i = 0; i < n.size(); ++i)
                        if(n[i] < 0)
                            state.setSpeciesAmount(i, options.learning.epsilon, "mol");

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

                    result.timing.estimate_database_priority_update = toc(PRIORITY_UPDATE_STEP);

                    // Mark the estimated state as accepted
                    result.estimate.accepted = true;

                    return;
                }
            }
        }

        result.estimate.accepted = false;
    }

    //=================================================================================================================
    //
    // MISCELLANEOUS METHODS
    //
    //=================================================================================================================

    /// Set the options of the smart equilibrium solver
    auto setOptions(SmartEquilibriumOptions const& opts) -> void
    {
        options = opts;
        solver.setOptions(opts.learning);
    }
};

SmartEquilibriumSolver::SmartEquilibriumSolver(ChemicalSystem const& system)
: pimpl(new Impl(EquilibriumSpecs::TP(system)))
{}

SmartEquilibriumSolver::SmartEquilibriumSolver(EquilibriumSpecs const& specs)
: pimpl(new Impl(specs))
{}

SmartEquilibriumSolver::SmartEquilibriumSolver(SmartEquilibriumSolver const& other)
: pimpl(new Impl(*other.pimpl))
{}

SmartEquilibriumSolver::~SmartEquilibriumSolver()
{}

auto SmartEquilibriumSolver::operator=(SmartEquilibriumSolver other) -> SmartEquilibriumSolver&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto SmartEquilibriumSolver::solve(ChemicalState& state) -> SmartEquilibriumResult
{
    return pimpl->solve(state);
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, EquilibriumRestrictions const& restrictions) -> SmartEquilibriumResult
{
    errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumRestrictions is currently not supported.");
    return {};
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, EquilibriumConditions const& conditions) -> SmartEquilibriumResult
{
    return pimpl->solve(state, conditions);
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions) -> SmartEquilibriumResult
{
    errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumRestrictions is currently not supported.");
    return {};
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, EquilibriumSensitivity& sensitivity) -> SmartEquilibriumResult
{
    errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumSensitivity is currently not supported.");
    return {};
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, EquilibriumRestrictions const& restrictions) -> SmartEquilibriumResult
{
    errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumSensitivity is currently not supported.");
    return {};
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, EquilibriumConditions const& conditions) -> SmartEquilibriumResult
{
    errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumSensitivity is currently not supported.");
    return {};
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions) -> SmartEquilibriumResult
{
    errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumSensitivity is currently not supported.");
    return {};
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, ArrayXdConstRef const& c0) -> SmartEquilibriumResult
{
    return pimpl->solve(state, c0);
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, EquilibriumRestrictions const& restrictions, ArrayXdConstRef const& c0) -> SmartEquilibriumResult
{
    errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumRestrictions is currently not supported.");
    return {};
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, EquilibriumConditions const& conditions, ArrayXdConstRef const& c0) -> SmartEquilibriumResult
{
    return pimpl->solve(state, conditions, c0);
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions, ArrayXdConstRef const& c0) -> SmartEquilibriumResult
{
    errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumRestrictions is currently not supported.");
    return {};
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, ArrayXdConstRef const& c0) -> SmartEquilibriumResult
{
    errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumSensitivity is currently not supported.");
    return {};
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, EquilibriumRestrictions const& restrictions, ArrayXdConstRef const& c0) -> SmartEquilibriumResult
{
    errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumSensitivity is currently not supported.");
    return {};
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, EquilibriumConditions const& conditions, ArrayXdConstRef const& c0) -> SmartEquilibriumResult
{
    errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumSensitivity is currently not supported.");
    return {};
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions, ArrayXdConstRef const& c0) -> SmartEquilibriumResult
{
    errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumSensitivity is currently not supported.");
    return {};
}

auto SmartEquilibriumSolver::setOptions(SmartEquilibriumOptions const& options) -> void
{
    pimpl->setOptions(options);
}

} // namespace Reaktoro
