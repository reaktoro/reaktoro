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

#include "KineticPathProfiler.hpp"

// C++ includes
#include <algorithm>
#include <iostream>

// Reaktoro includes
#include <Reaktoro/Common/OutputUtils.hpp>

namespace Reaktoro {

struct KineticPathProfiler::Impl
{
    /// The collected results of the smart kinetic path time step calculations.
    std::deque<SmartKineticResult> smart_results;

    /// The collected results of the kinetic path time step calculations.
    std::deque<KineticResult> results;

    /// The accumulated timing for the operations during kinetics calculations per time step.
    std::deque<SmartKineticTiming> timing_smart_kinetics_at_step;

    /// The accumulated timing for the operations during kinetics calculations per time step.
    std::deque<KineticTiming> timing_kinetics_at_step;

    /// The accumulated timing for the operations during equilibrium calculations per time step.
    std::deque<EquilibriumTiming> timing_equilibrium_at_step;

    /// The accumulated timing for the operations during smart equilibrium calculations per time step.
    std::deque<SmartEquilibriumTiming> timing_smart_equilibrium_at_step;

    /// The total accumulated timing for kinetics calculations.
    KineticTiming accumulated_timing_kinetics;

    /// The total accumulated timing for smart kinetics calculations.
    SmartKineticTiming accumulated_timing_smart_kinetics;

    /// The total accumulated timing for equilibrium calculations.
    EquilibriumTiming accumulated_timing_equilibrium;

    /// The total accumulated timing for smart equilibrium calculations.
    SmartEquilibriumTiming accumulated_timing_smart_equilibrium;

    /// Construct an instance of KineticPathProfiler::Impl.
    Impl()
    {
    }

    /// Update the profiler with the result of the last kinetic path time step.
    auto update(const SmartKineticResult& result) -> void
    {
        // Add the result into the vector of results
        results.push_back({});
        smart_results.push_back(result);

        // Calculate the total times at one step
        timing_kinetics_at_step.push_back({});
        timing_smart_kinetics_at_step.push_back(result.timing);
        timing_equilibrium_at_step.push_back( result.equilibrium.timing );
        timing_smart_equilibrium_at_step.push_back( result.smart_equilibrium.timing );

        // Accumulate the total times at all steps
        accumulated_timing_kinetics += timing_kinetics_at_step.back();
        accumulated_timing_smart_kinetics += timing_smart_kinetics_at_step.back();
        accumulated_timing_equilibrium += timing_equilibrium_at_step.back();
        accumulated_timing_smart_equilibrium += timing_smart_equilibrium_at_step.back();
    }

    /// Update the profiler with the result of the last kinetic path time step.
    auto update(const KineticResult& result) -> void
    {
        // Add the result into the vector of results
        results.push_back(result);
        smart_results.push_back({});

        // Calculate the total times at one step
        timing_kinetics_at_step.push_back( result.timing );
        timing_smart_kinetics_at_step.push_back( {} );
        timing_equilibrium_at_step.push_back( result.equilibrium.timing );
        timing_smart_equilibrium_at_step.push_back( result.smart_equilibrium.timing );

        // Accumulate the total times at all steps
        accumulated_timing_kinetics += timing_kinetics_at_step.back();
        accumulated_timing_smart_kinetics += timing_smart_kinetics_at_step.back();
        accumulated_timing_equilibrium += timing_equilibrium_at_step.back();
        accumulated_timing_smart_equilibrium += timing_smart_equilibrium_at_step.back();
    }

    /// Return the computing costs of all operations during a reactive transport calculation.
    auto computingCostsPerTimeStep() const -> KineticPathAnalysis::ComputingCostsPerTimeStep
    {
        KineticPathAnalysis::ComputingCostsPerTimeStep info;

        const auto num_time_steps = results.size();

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

        for(Index i = 0; i < num_time_steps; ++i)
        {
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
        }

        return info;
    }

    /// Return a summary of the performance analysis of all kinetics calculations.
    auto kineticsAnalysis() const -> KineticPathAnalysis::KineticAnalysis
    {
        KineticPathAnalysis::KineticAnalysis info;
        info.timing = accumulated_timing_kinetics;
        return info;
    }

    /// Return a summary of the performance analysis of all smart kinetics calculations.
    auto smartKineticsAnalysis() const -> KineticPathAnalysis::SmartKineticAnalysis
    {
        KineticPathAnalysis::SmartKineticAnalysis info;
        info.timing = accumulated_timing_smart_kinetics;

        // Count the accepted smart kinetic estimates and required learning operations
        for(const auto& result : smart_results)
            if(result.estimate.accepted) ++info.num_smart_kinetic_accepted_estimates;
            else ++info.num_smart_kinetic_required_learnings;

        // Set the total number of kinetic calculations
        info.num_kinetic_calculations =
                info.num_smart_kinetic_accepted_estimates + info.num_smart_kinetic_required_learnings;

        // Set the success rate at which smart equilibrium estimates were accepted
        info.smart_kinetic_estimate_acceptance_rate =
                static_cast<double>(info.num_smart_kinetic_accepted_estimates) / info.num_kinetic_calculations;

        // The number of time steps (= the number of collected KineticPathResult objects)
        const auto num_time_steps = smart_results.size();

        // For each time step, identify where learning was required
        for(Index i = 0; i < num_time_steps; ++i)
        {
            // If chemical state wasn't accepted and learning had to be performed, add the step-index to the dequeue
            if(!smart_results[i].estimate.accepted)
                info.steps_where_learning_was_required.push_back(i);
        }

        return info;
    }

    /// Return a summary of the performance analysis of all equilibrium calculations.
    auto equilibriumAnalysis() const -> KineticPathAnalysis::EquilibriumAnalysis
    {
        KineticPathAnalysis::EquilibriumAnalysis info;
        info.timing = accumulated_timing_equilibrium;
        return info;
    }

    /// Return a summary of the performance analysis of all smart equilibrium calculations.
    auto smartEquilibriumAnalysis() const -> KineticPathAnalysis::SmartEquilibriumAnalysis
    {
        KineticPathAnalysis::SmartEquilibriumAnalysis info;
        info.timing = accumulated_timing_smart_equilibrium;

        // Count the accepted smart equilibrium estimates and required learning operations
        for(const auto& result : results)
            if(result.smart_equilibrium.estimate.accepted) ++info.num_smart_equilibrium_accepted_estimates;
            else ++info.num_smart_equilibrium_required_learnings;

        // Set the total number of equilibrium calculations
        info.num_equilibrium_calculations =
            info.num_smart_equilibrium_accepted_estimates + info.num_smart_equilibrium_required_learnings;

        // Set the success rate at which smart equilibrium estimates were accepted
        info.smart_equilibrium_estimate_acceptance_rate =
            static_cast<double>(info.num_smart_equilibrium_accepted_estimates) / info.num_equilibrium_calculations;

        // The number of time steps (= the number of collected KineticPathResult objects)
        const auto num_time_steps = results.size();

        // For each time step, identify where learning was required
        for(Index i = 0; i < num_time_steps; ++i)
        {
            // If chemical state wasn't accepted and learning had to be performed, add the step-index to the dequeue
            if(!results[i].smart_equilibrium.estimate.accepted)
                    info.steps_where_learning_was_required.push_back(i);
        }

        return info;
    }

    /// Return the performance analysis of all operations in the reactive transport simulation.
    auto analysis() const -> KineticPathAnalysis
    {
        KineticPathAnalysis info;
        info.computing_costs_per_time_step = computingCostsPerTimeStep();
        info.kinetics = kineticsAnalysis();
        info.smart_kinetics = smartKineticsAnalysis();
        info.equilibrium = equilibriumAnalysis();
        info.smart_equilibrium = smartEquilibriumAnalysis();
        return info;
    }

};

KineticPathProfiler::KineticPathProfiler()
: pimpl(new Impl())
{}

KineticPathProfiler::KineticPathProfiler(const KineticPathProfiler& other)
: pimpl(new Impl(*other.pimpl))
{}

KineticPathProfiler::~KineticPathProfiler()
{}

auto KineticPathProfiler::operator=(KineticPathProfiler other) -> KineticPathProfiler&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto KineticPathProfiler::update(const KineticResult& result) -> void
{
    pimpl->update(result);
}

auto KineticPathProfiler::update(const SmartKineticResult& result) -> void
{
    pimpl->update(result);
}

auto KineticPathProfiler::analysis() const -> KineticPathAnalysis
{
    return pimpl->analysis();
}

auto KineticPathProfiler::results() const -> const std::deque<KineticResult>&
{
    return pimpl->results;
}

auto operator<<(std::ostream& out, const KineticPathProfiler& prof) -> std::ostream&
{
    auto seconds = [](double x) { return std::to_string(x) + "s"; };
    auto percent = [](double x, double y, std::string msg) { return std::to_string(x/y * 100) + "% " + msg; };
    auto status = [&](double x, double y, std::string msg) { return percent(x, y, msg) + " (" + seconds(x) + ")"; };

    const auto& analysis = prof.analysis();
    const auto& kin_timing = analysis.smart_kinetics.timing;

    out << "# ----------------------------------------------------------------------------------" << std::endl;
    out << "# Computing costs analysis of the operations in smart chemical kinetic calculations" << std::endl;
    out << "# ----------------------------------------------------------------------------------" << std::endl;
    out << "# solve                         = " << seconds(analysis.smart_equilibrium.timing.solve) << std::endl;
    out << "#   learning                      = " << status(kin_timing.learn, kin_timing.solve, "of total solve time") << std::endl;
    out << "#     integration                   = " << status(kin_timing.learn_integration, kin_timing.learn, "of learning time") << std::endl;
    out << "#       equilibration                 = " << status(kin_timing.learn_equilibration, kin_timing.learn, "of learning time") << std::endl;
    out << "#       chemical_properties           = " << status(kin_timing.learn_chemical_properties, kin_timing.learn, "of learning time") << std::endl;
    out << "#       reaction rates                = " << status(kin_timing.learn_reaction_rates, kin_timing.learn, "of learning time") << std::endl;
    out << "#     error control matrices        = " << status(kin_timing.learn_error_control_matrices, kin_timing.learn, "of learning time") << std::endl;
    out << "#     storage                       = " << status(kin_timing.learn_storage, kin_timing.learn, "of learning time") << std::endl;
    out << "#   estimate                      = " << status(kin_timing.estimate, kin_timing.solve, "of total solve time") << std::endl;
    out << "#     search                        = " << status(kin_timing.estimate_search, kin_timing.estimate, "of estimate time") << std::endl;
    out << "#       error control                 = " << status(kin_timing.estimate_error_control, kin_timing.estimate, "of estimate time") << std::endl;
    out << "#     taylor                        = " << status(kin_timing.estimate_taylor, kin_timing.estimate, "of estimate time") << std::endl;
    out << "#     database update               = " << status(kin_timing.estimate_database_priority_update, kin_timing.estimate, "of estimate time") << std::endl;
    out << "# ----------------------------------------------------------------------------------" << std::endl;
    out << "#" << std::endl;
    out << "# ----------------------------------------------------------------------------------" << std::endl;
    out << "# Overall computing costs in all smart chemical kinetic calculations" << std::endl;
    out << "# ----------------------------------------------------------------------------------" << std::endl;
    out << "# number of kinetics calculations             = " << analysis.smart_kinetics.num_kinetic_calculations << std::endl;
    out << "# number of smart kinetics accepted estimates = " << analysis.smart_kinetics.num_smart_kinetic_accepted_estimates << std::endl;
    out << "# number of smart kinetics required learnings = " << analysis.smart_kinetics.num_smart_kinetic_required_learnings << std::endl;
    out << "# smart kinetics estimate acceptance rate     = " << analysis.smart_kinetics.smart_kinetic_estimate_acceptance_rate * 100 << "%" << std::endl;
    out << "# ----------------------------------------------------------------------------------" << std::endl;
//    out << "#" << std::endl;
//    out << "# ----------------------------------------------------------------------------------" << std::endl;
//    out << "# Number of time steps number where learning was required      " << std::endl;
//    out << "# ----------------------------------------------------------------------------------" << std::endl;
//    const auto num_time_steps = analysis.smart_kinetics.steps_where_learning_was_required.size();
//    for(Index i = 0; i < num_time_steps; ++i)
//        out << analysis.smart_kinetics.steps_where_learning_was_required[i] << ", ";
//    out << std::endl;
//    out << "# ----------------------------------------------------------------------------------" << std::endl;
//    out << std::endl;

    const auto& timing = analysis.smart_equilibrium.timing;


//    out << "# -------------------------------------------------------------------------------------" << std::endl;
//    out << "# Computing costs analysis of the operations in smart chemical equilibrium calculations" << std::endl;
//    out << "# -------------------------------------------------------------------------------------" << std::endl;
//    out << "# solve                         = " << seconds(analysis.smart_equilibrium.timing.solve) << std::endl;
//    out << "#   learning                      = " << status(timing.learn, timing.solve, "of total solve time") << std::endl;
//    out << "#     gibbs_energy_minimization     = " << status(timing.learn_gibbs_energy_minimization, timing.learn, "of learning time") << std::endl;
//    out << "#     chemical_properties           = " << status(timing.learn_chemical_properties, timing.learn, "of learning time") << std::endl;
//    out << "#     sensitivity_matrix            = " << status(timing.learn_sensitivity_matrix, timing.learn, "of learning time") << std::endl;
//    out << "#     storage                       = " << status(timing.learn_storage, timing.learn, "of learning time") << std::endl;
//    out << "#   estimate                      = " << status(timing.estimate, timing.solve, "of total solve time") << std::endl;
//    out << "#     search                        = " << status(timing.estimate_search, timing.estimate, "of estimate time") << std::endl;
//    out << "#     taylor                        = " << status(timing.estimate_taylor, timing.estimate, "of estimate time") << std::endl;
//    out << "#     error control                 = " << status(timing.estimate_error_control, timing.estimate, "of estimate time") << std::endl;
//    out << "# -------------------------------------------------------------------------------------" << std::endl;
//    out << "#" << std::endl;
//    out << "# -------------------------------------------------------------------------------------" << std::endl;
//    out << "# Overall computing costs in all smart chemical equilibrium calculations" << std::endl;
//    out << "# -------------------------------------------------------------------------------------" << std::endl;
//    out << "# number of equilibrium calculations             = " << analysis.smart_equilibrium.num_equilibrium_calculations << std::endl;
//    out << "# number of smart equilibrium accepted estimates = " << analysis.smart_equilibrium.num_smart_equilibrium_accepted_estimates << std::endl;
//    out << "# number of smart equilibrium required learnings = " << analysis.smart_equilibrium.num_smart_equilibrium_required_learnings << std::endl;
//    out << "# smart equilibrium estimate acceptance rate     = " << analysis.smart_equilibrium.smart_equilibrium_estimate_acceptance_rate * 100 << "%" << std::endl;
//    out << "# -------------------------------------------------------------------------------------" << std::endl;
//    out << "#" << std::endl;
//    out << "# -------------------------------------------------------------------------------------" << std::endl;
//    out << "# Number of time steps number where learning was required      " << std::endl;
//    out << "#-------------------------------------------------------------------------------------" << std::endl;
//    for(Index i = 0; i < num_time_steps; ++i)
//        out << analysis.smart_equilibrium.steps_where_learning_was_required[i] << ", ";
//    out << std::endl;
//    out << "# -------------------------------------------------------------------------------------" << std::endl;
//    out << std::endl;


    return out;
}

} // namespace Reaktoro
