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
#include <Reaktoro/Common/OutputUtils.hpp>

namespace Reaktoro {
namespace internal {

template<typename ResultType>
auto accumulateTimingsFromResults(const std::vector<ResultType>& results)
{
    using TimingType = decltype(results.front().timing);
    TimingType accumulated_timing = TimingType{};
    if(results.empty()) return accumulated_timing;
    for(const auto& result : results)
        accumulated_timing += result.timing;
    return accumulated_timing;
}

} // namespace internal

struct ReactiveTransportProfiler::Impl
{
    /// The collected results of the reactive transport time step calculations.
    std::deque<ReactiveTransportResult> results;

    /// The accumulated timing for the operations during fluid element transport calculations per time step.
    std::deque<TransportTiming> timing_transport_at_step;

    /// The accumulated timing for the operations during equilibrium calculations per time step.
    std::deque<EquilibriumTiming> timing_equilibrium_at_step;

    /// The accumulated timing for the operations during smart equilibrium calculations per time step.
    std::deque<SmartEquilibriumTiming> timing_smart_equilibrium_at_step;

    /// The accumulated timing for the operations during kinetics calculations per time step.
    std::deque<KineticTiming> timing_kinetics_at_step;

    /// The accumulated timing for the operations during smart kinetics calculations per time step.
    std::deque<SmartKineticTiming> timing_smart_kinetics_at_step;

    /// The total accumulated timing for fluid element transport calculations.
    TransportTiming accumulated_timing_transport;

    /// The total accumulated timing for equilibrium calculations.
    EquilibriumTiming accumulated_timing_equilibrium;

    /// The total accumulated timing for smart equilibrium calculations.
    SmartEquilibriumTiming accumulated_timing_smart_equilibrium;

    /// The total accumulated timing for smart equilibrium calculations.
    KineticTiming accumulated_timing_kinetics;

    /// The total accumulated timing for smart kinetics calculations.
    SmartKineticTiming accumulated_timing_smart_kinetics;

    /// Construct an instance of ReactiveTransportProfiler::Impl.
    Impl()
    {
    }

    /// Update the profiler with the result of the last reactive transport time step.
    auto update(const ReactiveTransportResult& result) -> void
    {
        // Add the result into the vector of results
        results.push_back(result);

        // Calculate the total times (at one step) by summing the times in instances (elements or cells)
        timing_transport_at_step.push_back( internal::accumulateTimingsFromResults(result.transport_of_element) );
        timing_equilibrium_at_step.push_back( internal::accumulateTimingsFromResults(result.equilibrium_at_cell) );
        timing_smart_equilibrium_at_step.push_back( internal::accumulateTimingsFromResults(result.smart_equilibrium_at_cell) );
        timing_kinetics_at_step.push_back( internal::accumulateTimingsFromResults(result.kinetics_at_cell) );
        timing_smart_kinetics_at_step.push_back( internal::accumulateTimingsFromResults(result.smart_kinetics_at_cell) );

        // Accumulate the total times (at steps) to track the time over all cells
        accumulated_timing_transport += timing_transport_at_step.back();
        accumulated_timing_equilibrium += timing_equilibrium_at_step.back();
        accumulated_timing_smart_equilibrium += timing_smart_equilibrium_at_step.back();
        accumulated_timing_kinetics += timing_kinetics_at_step.back();
        accumulated_timing_smart_kinetics += timing_smart_kinetics_at_step.back();
    }

    /// Return the computing costs of all operations during a reactive transport calculation.
    auto computingCostsPerTimeStep() const -> ReactiveTransportAnalysis::ComputingCostsPerTimeStep
    {
        ReactiveTransportAnalysis::ComputingCostsPerTimeStep info;

        const auto num_time_steps = results.size();

        info.transport.resize(num_time_steps);
        info.equilibrium.resize(num_time_steps);

        info.smart_equilibrium.resize(num_time_steps);
        info.smart_equilibrium_with_ideal_search.resize(num_time_steps);

        info.smart_equilibrium_estimate.resize(num_time_steps);
        info.smart_equilibrium_search.resize(num_time_steps);
        info.smart_equilibrium_error_control.resize(num_time_steps);
        info.smart_equilibrium_taylor.resize(num_time_steps);
        info.smart_equilibrium_database_priority_update.resize(num_time_steps);

        info.smart_equilibrium_learn.resize(num_time_steps);
        info.smart_equilibrium_gibbs_energy_minimization.resize(num_time_steps);
        info.smart_equilibrium_chemical_properties.resize(num_time_steps);
        info.smart_equilibrium_sensitivity_matrix.resize(num_time_steps);
        info.smart_equilibrium_error_control_matrices.resize(num_time_steps);
        info.smart_equilibrium_storage.resize(num_time_steps);

        info.kinetics.resize(num_time_steps);
        info.kinetics_equilibration.resize(num_time_steps);
        info.kinetics_properties.resize(num_time_steps);
        info.kinetics_with_ideal_properties.resize(num_time_steps);

        info.smart_kinetics.resize(num_time_steps);
        info.smart_kinetics_with_ideal_search.resize(num_time_steps);

        info.smart_kinetics_estimate.resize(num_time_steps);
        info.smart_kinetics_search.resize(num_time_steps);
        info.smart_kinetics_error_control.resize(num_time_steps);
        info.smart_kinetics_taylor.resize(num_time_steps);

        info.smart_kinetics_learn.resize(num_time_steps);
        info.smart_kinetics_chemical_properties.resize(num_time_steps);
        info.smart_kinetics_equilibration.resize(num_time_steps);

        for(Index i = 0; i < num_time_steps; ++i)
        {
            info.transport[i] = timing_transport_at_step[i].step;
            info.equilibrium[i] = timing_equilibrium_at_step[i].solve;

            info.smart_equilibrium[i] = timing_smart_equilibrium_at_step[i].solve;
            info.smart_equilibrium_with_ideal_search[i] = info.smart_equilibrium[i] - timing_smart_equilibrium_at_step[i].estimate_search - timing_smart_equilibrium_at_step[i].estimate_database_priority_update;

            info.smart_equilibrium_estimate[i] = timing_smart_equilibrium_at_step[i].estimate;
            info.smart_equilibrium_search[i] = timing_smart_equilibrium_at_step[i].estimate_search;
            info.smart_equilibrium_error_control[i] = timing_smart_equilibrium_at_step[i].estimate_error_control;
            info.smart_equilibrium_taylor[i] = timing_smart_equilibrium_at_step[i].estimate_taylor;
            info.smart_equilibrium_database_priority_update[i] = timing_smart_equilibrium_at_step[i].estimate_database_priority_update;

            info.smart_equilibrium_learn[i] = timing_smart_equilibrium_at_step[i].learn;
            info.smart_equilibrium_gibbs_energy_minimization[i] = timing_smart_equilibrium_at_step[i].learn_gibbs_energy_minimization;
            info.smart_equilibrium_chemical_properties[i] = timing_smart_equilibrium_at_step[i].learn_chemical_properties;
            info.smart_equilibrium_sensitivity_matrix[i] = timing_smart_equilibrium_at_step[i].learn_sensitivity_matrix;
            info.smart_equilibrium_error_control_matrices[i] = timing_smart_equilibrium_at_step[i].learn_error_control_matrices;
            info.smart_equilibrium_storage[i] = timing_smart_equilibrium_at_step[i].learn_storage;

            info.kinetics[i] = timing_kinetics_at_step[i].solve;
            info.kinetics_equilibration[i] = timing_kinetics_at_step[i].integrate_equilibration;
            info.kinetics_properties[i] = timing_kinetics_at_step[i].integrate_chemical_properties;

            info.smart_kinetics[i] = timing_smart_kinetics_at_step[i].solve;
            info.smart_kinetics_with_ideal_search[i] = info.smart_kinetics[i] - timing_smart_kinetics_at_step[i].estimate_search;

            info.smart_kinetics_estimate[i] = timing_smart_kinetics_at_step[i].estimate;
            info.smart_kinetics_search[i] = timing_smart_kinetics_at_step[i].estimate_search;
            info.smart_kinetics_error_control[i] = timing_smart_kinetics_at_step[i].estimate_error_control;
            info.smart_kinetics_taylor[i] = timing_smart_kinetics_at_step[i].estimate_taylor;

            info.smart_kinetics_learn[i] = timing_smart_kinetics_at_step[i].learn;
            info.smart_kinetics_chemical_properties[i] = timing_smart_kinetics_at_step[i].learn_chemical_properties;
            info.smart_kinetics_equilibration[i] = timing_smart_kinetics_at_step[i].learn_equilibration;
        }

        return info;
    }

    /// Return a summary of the performance analysis of all transport calculations.
    auto transportAnalysis() const -> ReactiveTransportAnalysis::TransportAnalysis
    {
        ReactiveTransportAnalysis::TransportAnalysis info;
        info.timing = accumulated_timing_transport;
        return info;
    }

    /// Return a summary of the performance analysis of all equilibrium calculations.
    auto equilibriumAnalysis() const -> ReactiveTransportAnalysis::EquilibriumAnalysis
    {
        ReactiveTransportAnalysis::EquilibriumAnalysis info;
        info.timing = accumulated_timing_equilibrium;
        return info;
    }

    /// Return a summary of the performance analysis of all equilibrium calculations.
    auto kineticsAnalysis() const -> ReactiveTransportAnalysis::KineticsAnalysis
    {
        ReactiveTransportAnalysis::KineticsAnalysis info;
        info.timing = accumulated_timing_kinetics;
        return info;
    }

    /// Return a summary of the performance analysis of all smart equilibrium calculations.
    auto smartEquilibriumAnalysis() const -> ReactiveTransportAnalysis::SmartEquilibriumAnalysis
    {
        ReactiveTransportAnalysis::SmartEquilibriumAnalysis info;
        info.timing = accumulated_timing_smart_equilibrium;

        // Count the accepted smart equilibrium estimates and required learning operations
        for(const auto& result : results)
            for(const auto& smart_equilibrium_result : result.smart_equilibrium_at_cell)
                if(smart_equilibrium_result.estimate.accepted) ++info.num_smart_equilibrium_accepted_estimates;
                else ++info.num_smart_equilibrium_required_learnings;

        // Set the total number of equilibrium calculations
        info.num_equilibrium_calculations =
            info.num_smart_equilibrium_accepted_estimates + info.num_smart_equilibrium_required_learnings;

        // Set the success rate at which smart equilibrium estimates were accepted
        info.smart_equilibrium_estimate_acceptance_rate =
            static_cast<double>(info.num_smart_equilibrium_accepted_estimates) / info.num_equilibrium_calculations;

        // The number of time steps (= the number of collected ReactiveTransportResult objects)
        const auto num_time_steps = results.size();

        // Set the table of indicators that shows where the smart equilibrium estimate was accepted at a cell in a time step.
        info.cells_where_learning_was_required_at_step.resize(num_time_steps);

        // For each time step, identify the cells where learning was required
        for(Index i = 0; i < num_time_steps; ++i)
        {
            // For each cell, check if the smart estimation was accepted at current time step number
            const auto num_cells = results[i].smart_equilibrium_at_cell.size();
            for(Index j = 0; j < num_cells; ++j)
                if(!results[i].smart_equilibrium_at_cell[j].estimate.accepted)
                    info.cells_where_learning_was_required_at_step[i].push_back(j);
        }

        return info;
    }

    /// Return a summary of the performance analysis of all smart equilibrium calculations.
    auto smartKineticsAnalysis() const -> ReactiveTransportAnalysis::SmartKineticsAnalysis
    {
        ReactiveTransportAnalysis::SmartKineticsAnalysis info;
        info.timing = accumulated_timing_smart_kinetics;

        // Count the accepted smart equilibrium estimates and required learning operations
        for(const auto& result : results)
            for(const auto& smart_kinetics_result : result.smart_kinetics_at_cell)
                if(smart_kinetics_result.estimate.accepted) ++info.num_smart_kinetics_accepted_estimates;
                else ++info.num_smart_kinetics_required_learnings;

        // Set the total number of equilibrium calculations
        info.num_kinetics_calculations =
                info.num_smart_kinetics_accepted_estimates + info.num_smart_kinetics_required_learnings;

        // Set the success rate at which smart equilibrium estimates were accepted
        info.smart_kinetics_estimate_acceptance_rate =
                static_cast<double>(info.num_smart_kinetics_accepted_estimates) / info.num_kinetics_calculations;


        // The number of time steps (= the number of collected ReactiveTransportResult objects)
        const auto num_time_steps = results.size();

        // Set the table of indicators that shows where the smart kinetics estimate was accepted at a cell in a time step.
        info.cells_where_learning_was_required_at_step.resize(num_time_steps);

        // For each time step, identify the cells where learning was required
        for(Index i = 0; i < num_time_steps; ++i)
        {
            // For each cell, check if the smart estimation was accepted at current time step number
            const auto num_cells = results[i].smart_kinetics_at_cell.size();
            for(Index j = 0; j < num_cells; ++j)
                if (!results[i].smart_kinetics_at_cell[j].estimate.accepted)
                    info.cells_where_learning_was_required_at_step[i].push_back(j);

        }

        return info;
    }

    /// Return the performance analysis of all operations in the reactive transport simulation.
    auto analysis() const -> ReactiveTransportAnalysis
    {
        ReactiveTransportAnalysis info;
        info.computing_costs_per_time_step = computingCostsPerTimeStep();
        info.transport = transportAnalysis();
        info.equilibrium = equilibriumAnalysis();
        info.smart_equilibrium = smartEquilibriumAnalysis();
        info.kinetics = kineticsAnalysis();
        info.smart_kinetics = smartKineticsAnalysis();

        return info;
    }
};

ReactiveTransportProfiler::ReactiveTransportProfiler()
: pimpl(new Impl())
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

auto ReactiveTransportProfiler::update(const ReactiveTransportResult& result) -> void
{
    pimpl->update(result);
}

auto ReactiveTransportProfiler::analysis() const -> ReactiveTransportAnalysis
{
    return pimpl->analysis();
}

auto ReactiveTransportProfiler::results() const -> const std::deque<ReactiveTransportResult>&
{
    return pimpl->results;
}

// auto operator<<(std::ostream& out, const ReactiveTransportProfiler::SmartEquilibriumProfiling& prof) -> std::ostream&
// {
//     const auto& timing = prof.timing;

//     auto seconds = [](double x) { return std::to_string(x) + "s"; };
//     auto percent = [](double x, double y, std::string msg) { return std::to_string(x/y * 100) + "% " + msg; };
//     auto status = [&](double x, double y, std::string msg) { return percent(x, y, msg) + " (" + seconds(x) + ")"; };

//     out << "# -------------------------------------------------------------------------------------" << std::endl;
//     out << "# Computing costs analysis of the operations in smart chemical equilibrium calculations" << std::endl;
//     out << "# -------------------------------------------------------------------------------------" << std::endl;
//     out << "# solve                         = " << seconds(timing.solve) << std::endl;
//     out << "#   learning                      = " << status(timing.learning, timing.solve, "of total solve time") << std::endl;
//     out << "#     gibbs_energy_minimization     = " << status(timing.learn_gibbs_energy_minimization, timing.learning, "of learning time") << std::endl;
//     out << "#     chemical_properties           = " << status(timing.learn_chemical_properties, timing.learning, "of learning time") << std::endl;
//     out << "#     sensitivity_matrix            = " << status(timing.learn_sensitivity_matrix, timing.learning, "of learning time") << std::endl;
//     out << "#     storage                       = " << status(timing.learning_storage, timing.learning, "of learning time") << std::endl;
//     out << "#   estimate                      = " << status(timing.estimate, timing.solve, "of total solve time") << std::endl;
//     out << "#     search                        = " << status(timing.estimate_search, timing.estimate, "of estimate time") << std::endl;
//     out << "#     mat_vec_mul                   = " << status(timing.estimate_taylor, timing.estimate, "of estimate time") << std::endl;
//     out << "#     acceptance                    = " << status(timing.estimate_error_control, timing.estimate, "of estimate time") << std::endl;
//     out << "# -------------------------------------------------------------------------------------" << std::endl;
//     out << "#" << std::endl;
//     out << "# ----------------------------------------------------------------------" << std::endl;
//     out << "# Overall computing costs in all smart chemical equilibrium calculations" << std::endl;
//     out << "# ----------------------------------------------------------------------" << std::endl;
//     out << "# number of equilibrium calculations             = " << prof.num_equilibrium_calculations << std::endl;
//     out << "# number of smart equilibrium accepted estimates = " << prof.num_smart_equilibrium_accepted_estimates << std::endl;
//     out << "# number of smart equilibrium required learnings = " << prof.num_smart_equilibrium_required_learnings << std::endl;
//     out << "# smart equilibrium estimate acceptance rate     = " << prof.smart_equilibrium_estimate_acceptance_rate * 100 << "%" << std::endl;
//     out << "# ----------------------------------------------------------------------" << std::endl;
//     out << "#" << std::endl;
//     out << "# -------------------------------------------------------------------------" << std::endl;
//     out << "# Each row below contains a time step number followed by cell indices where" << std::endl;
//     out << "# learning was required because the smart estimation was not accepted      " << std::endl;
//     out << "# -------------------------------------------------------------------------" << std::endl;
//     const auto num_time_steps = prof.cells_where_learning_was_required_at_step.size();
//     for(Index i = 0; i < num_time_steps; ++i)
//         out << i << "," << prof.cells_where_learning_was_required_at_step[i] << std::endl;
// }

// auto operator<<(std::ostream& out, const ReactiveTransportProfiler::ComputingCostsPerTimeStep& costs) -> std::ostream&
// {
//     // Print the names of the collected computing costs in the header
//     out << "timestep" << ",";
//     out << "transport" << ",";
//     out << "equilibrium" << ",";
//     out << "smart_equilibrium" << ",";
//     out << "smart_equilibrium_with_ideal_search" << ",";
//     out << "smart_equilibrium_estimate" << ",";
//     out << "smart_equilibrium_search" << ",";
//     out << "smart_equilibrium_gibbs_energy_minimization" << ",";
//     out << "smart_equilibrium_storage";
//     out << std::endl;

//     // Print the values of the collected computing costs in each subsequent row
//     for(std::size_t i = 0; i < costs.transport.size(); ++i)
//     {
//         out << i << ",";
//         out << costs.transport[i] << ",";
//         out << costs.equilibrium[i] << ",";
//         out << costs.smart_equilibrium[i] << ",";
//         out << costs.smart_equilibrium_with_ideal_search[i] << ",";
//         out << costs.smart_equilibrium_estimate[i] << ",";
//         out << costs.smart_equilibrium_search[i] << ",";
//         out << costs.smart_equilibrium_gibbs_energy_minimization[i] << ",";
//         out << costs.smart_equilibrium_storage[i];
//         out << std::endl;
//     }

//     return out;
// }

} // namespace Reaktoro
