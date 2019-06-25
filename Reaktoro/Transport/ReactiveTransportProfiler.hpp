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

#pragma once

// C++ includes
#include <string>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>

namespace Reaktoro {

// Forward declarations
struct Results;
struct ReactiveTransportResult;

/// Use this class to postprocess and analyse accumulated results of reactive transport.
class ReactiveTransportProfiler
{
public:
    /// Construct ReactiveTransportResult instance
    ReactiveTransportProfiler(
        const std::string& results_folder,
        const std::string& file,
        const bool& smart);

    /// Process results collected on one step of reactive transport
    /// \param rt_results that stores the results of the reactive transport
    /// @see ReactiveTransportResult
    auto process(ReactiveTransportResult& rt_result) -> void;

    /// Update the file with results results collected on one step of reactive transport
    /// \param step number of step of reactive transport to be added to the file with results
    auto output(const Index& step) -> void;

    /// Summarize the profiling results of reactive transport
    auto summarize(Results& results) -> void;

private:
    /// Name of the file and folder with a status output
    std::string folder;
    std::string file;

    /// Used to store statistics information about the smart equilibrium algorithm.
    struct Statistics
    {
        /// Total time for search operations
        double time_estimate = 0.0;

        /// Time for search operations (part of estimation)
        double time_search = 0.0;

        /// Time for matrix-vector multiplications (part of estimation)
        double time_mat_vect_mult = 0.0;

        /// Time for acceptance test (part of estimation)
        double time_acceptance = 0.0;

        /// Total time for learn operations
        double time_learn = 0.0;

        /// Time for store operations (part of learning)
        double time_store = 0.0;

        /// Time for search operations (part of learning)
        double time_gibbs_min = 0.0;

        /// The size of the search tree
        Index tree_size = 0;

        /// Counter for the smart statuses
        int learning_counter = 0;

        /// Counter for the smart statuses
        int total_counter = 0;

    };

    /// The vector of the  CPU times for learning and estimating time
    Statistics total_stats;

    /// The vector of the  CPU times for learning and estimating time
    std::vector<Statistics> step_stats;

    /// The vector of the  CPU times for learning and estimating time
    std::vector<bool> statuses;

    /// Vector of reactive transport times and equilibrium times
    std::vector<double> rt_times;
    std::vector<double> eq_times;

    /// Flag whether smart of conventional solver was used
    bool smart;
};

} // namespace Reaktoro