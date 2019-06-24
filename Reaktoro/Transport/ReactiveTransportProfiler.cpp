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
#include <fstream>
#include <iomanip>

// Eigen includes
#include <Reaktoro/deps/eigen3/Eigen/Dense>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Profiling.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Transport/ReactiveTransportSolver.hpp>

namespace Reaktoro {

/// Implementation of a ReactiveTranportProfiler that collects information accumulated during the reactive transport
ReactiveTransportProfiler::ReactiveTransportProfiler(
    const std::string& folder,
    const std::string& file,
    const bool& smart)
: folder(folder), file(file), smart(smart)
{}

/// Process results collected per each step of reactive transport
auto ReactiveTransportProfiler::process(ReactiveTransportResult& rt_result) -> void
{
    if(smart)
    {
        SmartEquilibriumResult smart_res = rt_result.equilibrium.smart;
        Statistics stats;

        // Initialize the results step-wise
        stats.time_estimate = smart_res.estimate_stats.time_estimate;
        stats.time_search = smart_res.estimate_stats.time_search;
        stats.time_mat_vect_mult = smart_res.estimate_stats.time_mat_vect_mult;
        stats.time_acceptance = smart_res.estimate_stats.time_acceptance;
        stats.time_learn = smart_res.learn_stats.time_learn;
        stats.time_store = smart_res.learn_stats.time_store;
        stats.time_gibbs_min = smart_res.learn_stats.time_gibbs_min;
        stats.tree_size = smart_res.tree_size;
        stats.learning_counter = (smart_res.learning_states_indx.empty()) ? 0 : smart_res.learning_states_indx.size();
        step_stats.emplace_back(stats);

        // Accumulate the step-wise results
        total_stats.time_estimate += smart_res.estimate_stats.time_estimate;
        total_stats.time_search += smart_res.estimate_stats.time_search;
        total_stats.time_mat_vect_mult += smart_res.estimate_stats.time_mat_vect_mult;
        total_stats.time_acceptance += smart_res.estimate_stats.time_acceptance;
        total_stats.time_learn += smart_res.learn_stats.time_learn;
        total_stats.time_store += smart_res.learn_stats.time_store;
        total_stats.time_gibbs_min += smart_res.learn_stats.time_gibbs_min;
        total_stats.tree_size += smart_res.tree_size;
        total_stats.learning_counter += stats.learning_counter ;
        total_stats.total_counter += rt_result.ncells;

        // Fill in the statuses array based on the
        auto it = smart_res.learning_states_indx.begin();

        for(Index i = 0; i < rt_result.ncells; i++)
        {
            if(it != smart_res.learning_states_indx.end() && *it == i)
            {
                statuses.emplace_back(false);
                it++;
            }
            else
                statuses.emplace_back(true);
        }
    }
    else
    {
        EquilibriumResult res = rt_result.equilibrium;
        Statistics stats;

        // Initialize the results step-wise
        stats.time_learn = res.stats.time_learn;
        step_stats.emplace_back(stats);

        // Accumulate the step-wise results
        total_stats.time_learn += res.stats.time_learn;
    }

    // Add reactive transport and equilibrium time
    rt_times.emplace_back(rt_result.rt_time);
    eq_times.emplace_back(rt_result.eq_time);

    // Null the CUP time
    rt_result = ReactiveTransportResult();
}

/// Output results collected at each step of reactive transport
auto ReactiveTransportProfiler::output(const Index &step) -> void
{
    // The output stream of the data file
    std::ofstream datafile;

    // The floating-point precision in the output.
    int precision = 6;

    // Statuses for opening of the file
    auto opt = std::ofstream::out;

    // Depending on the step, open new file or append to existing one
    if(step != 0 && !file.empty()) opt |= std::ofstream::app;
    else if(step == 0 && !file.empty()) opt |= std::ofstream::trunc;

    if(smart)
    {
        // Document statuses
        // ----------------------------------------------------------------------
        datafile.open(folder + "/statuses.txt", opt);
        // Output statuses such that steps' go vertically and cells horizontally
        if(datafile.is_open())
        {
            // Output statuses collected while stepping with ReactiveTransportSolver
            for(bool st : statuses) datafile << std::to_string(st) << "\t";
            datafile << "\n";
        }
        statuses.clear();

        // Close the data file
        datafile.close();
    }
}

/// Summarize the profiling results of reactive transport
auto ReactiveTransportProfiler::summarize(Results& results) -> void
{
    if(smart)
    {
        std::cout << std::setprecision(6);

        std::cout << "\nlearning time                  : " << total_stats.time_learn << " (" << total_stats.time_learn / (total_stats.time_learn + total_stats.time_estimate) * 100 << "% of total time)\n";
        std::cout << "   - gibbs minimization time   : " << total_stats.time_gibbs_min << " (" << total_stats.time_gibbs_min / total_stats.time_learn * 100 << "% of learning time)\n";
        std::cout << "   - store time                : " << total_stats.time_store << " (" << total_stats.time_store / total_stats.time_learn * 100 << "% of learning time)\n";
        std::cout << "estimating time                : " << total_stats.time_estimate  << " (" << total_stats.time_estimate / (total_stats.time_learn + total_stats.time_estimate) * 100 << "% of total time)\n";
        std::cout << "   - search time               : " << total_stats.time_search << " (" << total_stats.time_search / total_stats.time_estimate * 100 << "% of estimate time)\n";
        std::cout << "   - matrix-vector mult. time  : " << total_stats.time_mat_vect_mult << " (" << total_stats.time_mat_vect_mult / total_stats.time_estimate * 100 << "% of estimate time)\n";
        std::cout << "   - acceptance test time      : " << total_stats.time_acceptance << " (" << total_stats.time_acceptance / total_stats.time_estimate * 100 << "% of estimate time)\n\n";

        results.smart_total = total_stats.time_learn + total_stats.time_estimate;
        results.smart_total_ideal_search = total_stats.time_learn + total_stats.time_estimate - total_stats.time_search;
        results.smart_total_ideal_search_store = total_stats.time_learn + total_stats.time_estimate - total_stats.time_search - total_stats.time_store;

        std::cout << "total time (ideal search)               : " << total_stats.time_learn + total_stats.time_estimate - total_stats.time_search << "\n";
        std::cout << "total time (ideal store + ideal search) : " << total_stats.time_learn + total_stats.time_estimate - total_stats.time_search - total_stats.time_store << "\n\n";
        std::cout << "total time                              : " << results.smart_total << "\n\n";

        std::cout << std::setprecision(2)
              << 100.00 *
                 double(total_stats.learning_counter) /
                 double(total_stats.total_counter ) << "% of training : "
              << total_stats.learning_counter << " cells out of "
              << total_stats.total_counter << "\n\n";
    }
    else
    {
        results.conv_total = total_stats.time_learn + total_stats.time_estimate;

        std::cout << std::setprecision(6)
                  <<  "total time                     : "
                  << results.conv_total << "\n\n";
    }

    // The output stream of the data file
    std::ofstream datafile;

    // The floating-point precision in the output.
    int precision = 6;

    // Document times
    // ----------------------------------------------------------------------
    // ----------------------------------------------------------------------
    datafile.open(folder + "/" + file + "-steps.txt", std::ofstream::out | std::ofstream::trunc);
    // Output statuses such that steps' go vertically and cells horizontally
    if(datafile.is_open())
    {
        /*
        if(smart)
        {
            datafile
                    << "time_transport   time_equilib    time_estimate   time_search     "
                    << "time_mat_vect   time_acceptance time_learn      time_store      "
                    << "time_gibbs_min tree   smart_counter\n";
        }else{
            datafile
                    << "time_transport   time_equilib    time_learn\n";
        }
        */
        datafile << std::scientific << std::setprecision(precision);
        unsigned length = rt_times.size();
        for(unsigned i = 0; i < length; i++)
        {
            datafile << rt_times[i] << "\t ";
            datafile << eq_times[i] << "\t ";

            if(smart)
            {
                // Output statuses collected while stepping with ReactiveTransportSolver
                datafile << step_stats[i].time_estimate << "\t "
                          << step_stats[i].time_search << "\t "
                          << step_stats[i].time_mat_vect_mult << "\t "
                          << step_stats[i].time_acceptance << "\t "
                          << step_stats[i].time_learn << "\t "
                          << step_stats[i].time_store << "\t "
                          << step_stats[i].time_gibbs_min << "\t"
                          << step_stats[i].tree_size << "\t"
                          << step_stats[i].learning_counter << "\n";

            } else
                datafile << step_stats[i].time_learn << "\n";
        }
    }

    // Close the data file
    datafile.close();
}

} // namespace Reaktoro
