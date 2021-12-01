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

#pragma once

// C++ includes
#include <deque>
#include <ostream>
#include <memory>

// Reaktoro includes
#include <Reaktoro/Transport/ReactiveTransportAnalysis.hpp>
#include <Reaktoro/Transport/ReactiveTransportResult.hpp>

namespace Reaktoro {

/// Provide mechanisms for analysing the operations in a reactive transport simulation.
class ReactiveTransportProfiler
{
public:
    /// Construct a default instance of ReactiveTransportProfiler.
    ReactiveTransportProfiler();

    /// Construct a copy of a ReactiveTransportProfiler instance.
    ReactiveTransportProfiler(const ReactiveTransportProfiler& other);

    /// Destroy this ReactiveTransportProfiler instance.
    virtual ~ReactiveTransportProfiler();

    /// Assign a copy of an ReactiveTransportProfiler instance.
    auto operator=(ReactiveTransportProfiler other) -> ReactiveTransportProfiler&;

    /// Update the profiler with the result of the last reactive transport time step.
    auto update(const ReactiveTransportResult& result) -> void;

    /// Return the performance analysis of all operations in the reactive transport simulation.
    auto analysis() const -> ReactiveTransportAnalysis;

    /// Return all collected results of the reactive transport simulation.
    auto results() const -> const std::deque<ReactiveTransportResult>&;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
