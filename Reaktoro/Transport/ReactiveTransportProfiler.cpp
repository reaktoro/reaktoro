// This file is part of Reaktoro (https://reaktoro.org).
//
// Reaktoro is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// Reaktoro is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#include "ReactiveTransportProfiler.hpp"

// C++ includes
#include <algorithm>
#include <iomanip>
#include <iostream>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/OutputUtils.hpp>
#include <Reaktoro/Transport/ReactiveTransportSolver.hpp>

namespace Reaktoro {
namespace internal {

template<typename ResultType>
auto accumulateTimingsFromResults(const std::vector<ResultType>& results)
{
    using TimingType = decltype(results.front().timing);
    if(results.empty()) return TimingType{};
    TimingType accumulated_timing = results.front().timing;
    for(const auto& result : results)
        accumulated_timing += result.timing;
    return accumulated_timing;
}

} // namespace internal

struct ReactiveTransportProfiler::Impl
{
    /// The reference to the reactive transport solver.
    const ReactiveTransportSolver& solver;

    /// The collected results of the reactive transport time step calculations.
    std::deque<ReactiveTransportResult> results;

    /// The accumulated timing for the operations during fluid element transport calculations per time step.
    std::deque<TransportTiming> timing_transport_at_step;

    /// The accumulated timing for the operations during equilibrium calculations per time step.
    std::deque<EquilibriumTiming> timing_equilibrium_at_step;

    /// The accumulated timing for the operations during smart equilibrium calculations per time step.
    std::deque<SmartEquilibriumTiming> timing_smart_equilibrium_at_step;

    /// The accumulated timing for the operations during fluid element transport calculations.
    AccumulatedTimings accumulated_timings;

    /// Construct an instance of ReactiveTransportProfiler::Impl.
    Impl(const ReactiveTransportSolver& solver)
    : solver(solver)
    {
    }

    /// Update the profiler with a new reactive transport time step profiling data.
    auto update() -> void
    {
        results.push_back(solver.result());

        const auto& last_result = results.back();

        timing_transport_at_step.push_back( internal::accumulateTimingsFromResults(last_result.transport_of_element) );
        timing_equilibrium_at_step.push_back( internal::accumulateTimingsFromResults(last_result.equilibrium_at_cell) );
        timing_smart_equilibrium_at_step.push_back( internal::accumulateTimingsFromResults(last_result.smart_equilibrium_at_cell) );

        accumulated_timings.transport += timing_transport_at_step.back();
        accumulated_timings.equilibrium += timing_equilibrium_at_step.back();
        accumulated_timings.smart_equilibrium += timing_smart_equilibrium_at_step.back();
    }

    /// Return the computing costs of all operations during a reactive transport calculation.
    auto computingCostsPerTimeStep() const -> ComputingCostsPerTimeStep
    {
        const auto num_time_steps = results.size();
        const auto dt = solver.timeStep();

        ComputingCostsPerTimeStep costs;
        costs.t.resize(num_time_steps);
        costs.transport.resize(num_time_steps);
        costs.equilibrium.resize(num_time_steps);
        costs.smart_equilibrium.resize(num_time_steps);
        costs.smart_equilibrium_with_ideal_search.resize(num_time_steps);
        costs.smart_equilibrium_estimate.resize(num_time_steps);
        costs.smart_equilibrium_nearest_neighbor_search.resize(num_time_steps);
        costs.smart_equilibrium_gibbs_energy_minimization.resize(num_time_steps);
        costs.smart_equilibrium_storage.resize(num_time_steps);

        double t = 0.0;

        for(Index i = 0; i < num_time_steps; ++i)
        {
            costs.t[i] = t;
            costs.transport[i] = timing_transport_at_step[i].step;
            costs.equilibrium[i] = timing_equilibrium_at_step[i].solve;
            costs.smart_equilibrium[i] = timing_smart_equilibrium_at_step[i].solve;
            costs.smart_equilibrium_with_ideal_search[i] = costs.smart_equilibrium[i] - timing_smart_equilibrium_at_step[i].estimate_search;
            costs.smart_equilibrium_estimate[i] = timing_smart_equilibrium_at_step[i].estimate;
            costs.smart_equilibrium_nearest_neighbor_search[i] = timing_smart_equilibrium_at_step[i].estimate_search;
            costs.smart_equilibrium_gibbs_energy_minimization[i] = timing_smart_equilibrium_at_step[i].learning_gibbs_energy_minimization;
            costs.smart_equilibrium_storage[i] = timing_smart_equilibrium_at_step[i].learning_storage;
            t += dt;
        }

        return costs;
    }

    /// Return a summary of the performance analysis of the smart equilibrium operations in a reactive transport calculation.
    auto smartEquilibriumProfiling() const -> SmartEquilibriumProfiling
    {
        SmartEquilibriumProfiling prof;

        prof.timing = accumulated_timings.smart_equilibrium;

        // Count the accepted smart equilibrium estimates and required learning operations
        for(const auto& result : results)
            for(const auto& smart_equilibrium_result : result.smart_equilibrium_at_cell)
                if(smart_equilibrium_result.estimate.accepted) ++prof.num_smart_equilibrium_accepted_estimates;
                else ++prof.num_smart_equilibrium_required_learnings;

        // Set the total number of equilibrium calculations
        prof.num_equilibrium_calculations =
            prof.num_smart_equilibrium_accepted_estimates +
                prof.num_smart_equilibrium_required_learnings;

        // Set the success rate at which smart equilibrium estimates were accepted.
        prof.smart_equilibrium_estimate_acceptance_rate =
            static_cast<double>(prof.num_smart_equilibrium_accepted_estimates) / prof.num_equilibrium_calculations;

        // The number of time steps (= the number of collected ReactiveTransportResult objects)
        const auto num_time_steps = results.size();

        // Set the table of indicators that show where smart equilibrium estimate was accepted at a cell in a time step.
        prof.cells_where_learning_was_required_at_step.resize(num_time_steps);

        // For each time step, identify the cells where learning was required
        for(Index i = 0; i < num_time_steps; ++i)
        {
            // For each cell, check if the smart estimation was accepted at current time step number
            const auto num_cells = results[i].smart_equilibrium_at_cell.size();
            for(Index j = 0; j < num_cells; ++j)
                if(results[i].smart_equilibrium_at_cell[j].estimate.accepted == false)
                    prof.cells_where_learning_was_required_at_step[i].push_back(j);
        }

        return prof;
    }

    /// Output the complete analysis of the performance of reactive transport simulation.
    auto output(std::string filename) -> void
    {

    }
};

ReactiveTransportProfiler::ReactiveTransportProfiler(const ReactiveTransportSolver& solver)
: pimpl(new Impl(solver))
{}

ReactiveTransportProfiler::ReactiveTransportProfiler(const ReactiveTransportProfiler& other)
: pimpl(new Impl(*other.pimpl))
{}

ReactiveTransportProfiler::~ReactiveTransportProfiler()
{}

auto ReactiveTransportProfiler::operator=(ReactiveTransportProfiler other) -> ReactiveTransportProfiler&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto ReactiveTransportProfiler::update() -> void
{
    pimpl->update();
}

auto ReactiveTransportProfiler::results() const -> const std::deque<ReactiveTransportResult>&
{
    return pimpl->results;
}

auto ReactiveTransportProfiler::accumulatedTimings() const -> AccumulatedTimings
{
    return pimpl->accumulated_timings;
}

auto ReactiveTransportProfiler::computingCostsPerTimeStep() const -> ComputingCostsPerTimeStep
{
    return pimpl->computingCostsPerTimeStep();
}

auto ReactiveTransportProfiler::smartEquilibriumProfiling() const -> SmartEquilibriumProfiling
{
    return pimpl->smartEquilibriumProfiling();
}

auto ReactiveTransportProfiler::output(std::string filename) -> void
{
    pimpl->output(filename);
}

auto operator<<(std::ostream& out, const ReactiveTransportProfiler::SmartEquilibriumProfiling& prof) -> std::ostream&
{
    const auto& timing = prof.timing;

    auto seconds = [](double x) { return std::to_string(x) + "s"; };
    auto percent = [](double x, double y, std::string msg) { return std::to_string(x/y * 100) + "% " + msg; };
    auto status = [&](double x, double y, std::string msg) { return percent(x, y, msg) + " (" + seconds(x) + ")"; };

    out << "# -------------------------------------------------------------------------------------" << std::endl;
    out << "# Computing costs analysis of the operations in smart chemical equilibrium calculations" << std::endl;
    out << "# -------------------------------------------------------------------------------------" << std::endl;
    out << "# solve                         = " << seconds(timing.solve) << std::endl;
    out << "#   learning                      = " << status(timing.learning, timing.solve, "of total solve time") << std::endl;
    out << "#     gibbs_energy_minimization     = " << status(timing.learning_gibbs_energy_minimization, timing.learning, "of learning time") << std::endl;
    out << "#     chemical_properties           = " << status(timing.learning_chemical_properties, timing.learning, "of learning time") << std::endl;
    out << "#     sensitivity_matrix            = " << status(timing.learning_sensitivity_matrix, timing.learning, "of learning time") << std::endl;
    out << "#     storage                       = " << status(timing.learning_storage, timing.learning, "of learning time") << std::endl;
    out << "#   estimate                      = " << status(timing.estimate, timing.solve, "of total solve time") << std::endl;
    out << "#     search                        = " << status(timing.estimate_search, timing.estimate, "of estimate time") << std::endl;
    out << "#     mat_vec_mul                   = " << status(timing.estimate_mat_vec_mul, timing.estimate, "of estimate time") << std::endl;
    out << "#     acceptance                    = " << status(timing.estimate_acceptance, timing.estimate, "of estimate time") << std::endl;
    out << "# -------------------------------------------------------------------------------------" << std::endl;
    out << "#" << std::endl;
    out << "# ----------------------------------------------------------------------" << std::endl;
    out << "# Overall computing costs in all smart chemical equilibrium calculations" << std::endl;
    out << "# ----------------------------------------------------------------------" << std::endl;
    out << "# number of equilibrium calculations             = " << prof.num_equilibrium_calculations << std::endl;
    out << "# number of smart equilibrium accepted estimates = " << prof.num_smart_equilibrium_accepted_estimates << std::endl;
    out << "# number of smart equilibrium required learnings = " << prof.num_smart_equilibrium_required_learnings << std::endl;
    out << "# smart equilibrium estimate acceptance rate     = " << prof.smart_equilibrium_estimate_acceptance_rate * 100 << "%" << std::endl;
    out << "# ----------------------------------------------------------------------" << std::endl;
    out << "#" << std::endl;
    out << "# -------------------------------------------------------------------------" << std::endl;
    out << "# Each row below contains a time step number followed by cell indices where" << std::endl;
    out << "# learning was required because the smart estimation was not accepted      " << std::endl;
    out << "# -------------------------------------------------------------------------" << std::endl;
    const auto num_time_steps = prof.cells_where_learning_was_required_at_step.size();
    for(Index i = 0; i < num_time_steps; ++i)
        out << i << "," << prof.cells_where_learning_was_required_at_step[i] << std::endl;
}

auto operator<<(std::ostream& out, const ReactiveTransportProfiler::ComputingCostsPerTimeStep& costs) -> std::ostream&
{
    // Print the names of the collected computing costs in the header
    out << "t" << ",";
    out << "transport" << ",";
    out << "equilibrium" << ",";
    out << "smart_equilibrium" << ",";
    out << "smart_equilibrium_with_ideal_search" << ",";
    out << "smart_equilibrium_estimate" << ",";
    out << "smart_equilibrium_nearest_neighbor_search" << ",";
    out << "smart_equilibrium_gibbs_energy_minimization" << ",";
    out << "smart_equilibrium_storage";
    out << std::endl;

    // Print the values of the collected computing costs in each subsequent row
    for(std::size_t i = 0; i < costs.t.size(); ++i)
    {
        out << costs.t[i] << ",";
        out << costs.transport[i] << ",";
        out << costs.equilibrium[i] << ",";
        out << costs.smart_equilibrium[i] << ",";
        out << costs.smart_equilibrium_with_ideal_search[i] << ",";
        out << costs.smart_equilibrium_estimate[i] << ",";
        out << costs.smart_equilibrium_nearest_neighbor_search[i] << ",";
        out << costs.smart_equilibrium_gibbs_energy_minimization[i] << ",";
        out << costs.smart_equilibrium_storage[i];
        out << std::endl;
    }

    return out;
}

} // namespace Reaktoro
