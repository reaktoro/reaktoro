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

/// Compute the step-round value of a given number.
/// This method computes the step-round value of a given number for a given step length separating the set of step-rounded values.
/// For example, if `step = 5` and `7.5 <= num <= 12.5`, `sround(num, step) => 10`.
/// @param num The number for which its step-rounded value is sought
/// @param step The length between step-rounded values
auto sround(double num, double step) -> long
{
    return round(num / step) * step;
}

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

/// Return the hash number for given step-rounded temperature (in K), step-rounded pressure (in Pa) and a vector of indices of primary species.
auto hash(long iT, long iP, ArrayXlConstRef const& iprimary) -> std::size_t
{
    std::hash<long> hasher;
    std::size_t seed = 2 + iprimary.size();
    seed ^= hasher(iT) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    seed ^= hasher(iP) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    for(auto const& i : iprimary)
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
        tic(SOLVE_STEP)

        // Reset the result of the last smart equilibrium calculation
        result = {};

        // Check input variables have all been set
        errorifnot(conditions.inputValues().allFinite(), "Ensure all input variables have been set in the EquilibriumConditions object.");

        // Perform a smart prediction of the chemical state
        timeit( predict(state, conditions), result.timing.prediction= )

        // Perform a learning step if the smart prediction is not satisfactory
        if(!result.prediction.accepted)
            timeit( learn(state, conditions), result.timing.learning= )

        result.timing.solve = toc(SOLVE_STEP);

        return result;
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
    // LEARN AND PREDICT METHODS
    //
    //=================================================================================================================

    /// Perform a learning operation in which a full chemical equilibrium calculation is performed.
    auto learn(ChemicalState& state, EquilibriumConditions const& conditions) -> void
    {
        //---------------------------------------------------------------------
        // GIBBS ENERGY MINIMIZATION CALCULATION DURING THE LEARNING PROCESS
        //---------------------------------------------------------------------
        tic(EQUILIBRIUM_STEP)

        // Perform a full chemical equilibrium solve with sensitivity derivatives calculation
        result.learning.solve = solver.solve(state, sensitivity, conditions);

        result.timing.learning_solve = toc(EQUILIBRIUM_STEP);

        //---------------------------------------------------------------------
        // STORAGE STEP DURING THE LEARNING PROCESS
        //---------------------------------------------------------------------
        tic(STORAGE_STEP)

        // Create an equilibrium predictor object with computed equilibrium state and its sensitivities
        EquilibriumPredictor predictor(state, sensitivity);

        // Generate the hash number for the temperature, pressure and indices of primary species in the state
        const auto iT = detail::sround(state.temperature(), options.temperature_step);
        const auto iP = detail::sround(state.pressure(), options.pressure_step);
        const auto iprimary = state.equilibrium().indicesPrimarySpecies();
        const auto label = detail::hash(iT, iP, iprimary);

        // std::cout << iT << " " << iP/1e5 << " " << iprimary.transpose() << label << "\n";

        // Find the index of the cluster that has same primary species
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

        result.timing.learning_storage = toc(STORAGE_STEP);
    }

    /// Perform a prediction operation in which a chemical equilibrium state is predicted using a first-order Taylor approximation.
    auto predict(ChemicalState& state, EquilibriumConditions const& conditions) -> void
    {
        // Set the prediction status to false at the beginning
        result.prediction.accepted = false;

        // Skip prediction operation if no cluster exists yet
        if(database.clusters.empty())
            return;

        const auto wvals = conditions.inputValues();
        const auto cvals = conditions.initialComponentAmountsGetOrCompute(state);

        const auto w = wvals.cast<double>();
        const auto c = cvals.cast<double>();

        // Auxiliary vectors used in the lambda function below to avoid repeated memory allocation
        VectorXd dw;
        VectorXd dc;

        // The function that checks if a record in the database pass the error test.
        auto pass_error_test = [&](Record const& record) mutable -> bool
        {
            // The primary species at the reference chemical state
            auto const& iprimary0 = record.state.equilibrium().indicesPrimarySpecies();

            // The equilibrium predictor calculator at the reference state
            auto const& predictor0 = record.predictor;

            const auto w0 = record.state.equilibrium().w();
            const auto c0 = record.state.equilibrium().c();

            dw = w - w0;
            dc = c - c0;

            using std::abs;

            for(auto ispecies : iprimary0)
            {
                const auto mu0 = predictor0.speciesChemicalPotentialReference(ispecies);
                const auto mu1 = predictor0.speciesChemicalPotentialPredicted(ispecies, dw, dc);
                if(abs(mu1 - mu0) >= options.reltol*abs(mu0) + options.abstol)
                    return false;
            }

            return true;
        };

        // Generate the hash number for the temperature, pressure and indices of primary species in the state
        const auto iT = detail::sround(state.temperature(), options.temperature_step);
        const auto iP = detail::sround(state.pressure(), options.pressure_step);
        const auto iprimary = state.equilibrium().indicesPrimarySpecies();
        const auto label = detail::hash(iT, iP, iprimary);

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

            // Iterate over all records in current cluster (using the order based on the priorities)
            for(auto irecord : records_ordering)
            {
                auto const& record = records[irecord];

                //---------------------------------------------------------------------
                // ERROR CONTROL STEP DURING THE ESTIMATE PROCESS
                //---------------------------------------------------------------------
                tic(ERROR_CONTROL_STEP)

                // Check if the current record passes the error test
                const auto success = pass_error_test(record);

                result.timing.prediction_error_control += toc(ERROR_CONTROL_STEP);

                if(success)
                {
                    //---------------------------------------------------------------------
                    // TAYLOR PREDICTION STEP DURING THE ESTIMATE PROCESS
                    //---------------------------------------------------------------------
                    tic(TAYLOR_STEP)

                    auto const& predictor0 = record.predictor;

                    predictor0.predict(state, conditions);

                    result.timing.prediction_taylor = toc(TAYLOR_STEP);

                    // Check if all projected species amounts are positive or at least very small negative values
                    auto const& n = state.speciesAmounts();

                    const double nmin = n.minCoeff();
                    const double nsum = n.sum();

                    if(nmin <= options.reltol_negative_amounts * nsum)
                        continue; // continue searching for a another record that produces positive amounts only or tolerable negative values

                    result.timing.prediction_search = toc(SEARCH_STEP);

                    //---------------------------------------------------------------------
                    // After the search is finished successfully
                    //---------------------------------------------------------------------

                    // Assign small positive values to all negative amounts
                    for(auto i = 0; i < n.size(); ++i)
                        if(n[i] < 0.0)
                            state.setSpeciesAmount(i, options.learning.epsilon);

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

                    result.timing.prediction_database_priority_update = toc(PRIORITY_UPDATE_STEP);

                    // Mark the predicted state as accepted
                    result.prediction.accepted = true;

                    return;
                }
            }
        }

        result.prediction.accepted = false;
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

auto SmartEquilibriumSolver::setOptions(SmartEquilibriumOptions const& options) -> void
{
    pimpl->setOptions(options);
}

} // namespace Reaktoro
