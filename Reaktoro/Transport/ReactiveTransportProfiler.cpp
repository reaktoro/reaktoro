// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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

#include "ReactiveTransportProfiler.hpp"

// C++ includes
#include <algorithm>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>

namespace Reaktoro {
namespace internal {

template<typename ResultType>
auto accumulateTimingsFromResults(const std::vector<ResultType>& results)
{
    using TimingType = decltype(results.front().timing);
    TimingType accumulated_timing = {};
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

    /// The total accumulated timing for fluid element transport calculations.
    TransportTiming accumulated_timing_transport;

    /// The total accumulated timing for equilibrium calculations.
    EquilibriumTiming accumulated_timing_equilibrium;

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

        // Accumulate the total times (at steps) to track the time over all cells
        accumulated_timing_transport += timing_transport_at_step.back();
        accumulated_timing_equilibrium += timing_equilibrium_at_step.back();
    }

    /// Return the computing costs of all operations during a reactive transport calculation.
    auto computingCostsPerTimeStep() const -> ReactiveTransportAnalysis::ComputingCostsPerTimeStep
    {
        ReactiveTransportAnalysis::ComputingCostsPerTimeStep info;

        const auto num_time_steps = results.size();

        info.transport.resize(num_time_steps);
        info.equilibrium.resize(num_time_steps);

        for(Index i = 0; i < num_time_steps; ++i)
        {
            info.transport[i] = timing_transport_at_step[i].step;
            info.equilibrium[i] = timing_equilibrium_at_step[i].solve;
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

    /// Return the performance analysis of all operations in the reactive transport simulation.
    auto analysis() const -> ReactiveTransportAnalysis
    {
        ReactiveTransportAnalysis info;
        info.computing_costs_per_time_step = computingCostsPerTimeStep();
        info.transport = transportAnalysis();
        info.equilibrium = equilibriumAnalysis();

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

} // namespace Reaktoro
