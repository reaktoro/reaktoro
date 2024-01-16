// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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
#include <chrono>

namespace Reaktoro {

using Time = std::chrono::time_point<std::chrono::high_resolution_clock>;

using Duration = std::chrono::duration<double>;

/// Return the time point now
/// @see elapsed
auto time() -> Time;

/// Return the elapsed time between two time points (in units of s)
/// @param end The end time point
/// @param end The begin time point
/// @return The elapsed time between *end* and *begin* in seconds
auto elapsed(Time const& end, Time const& begin) -> double;

/// Return the elapsed time between a time point and now (in units of s)
/// @param end The begin time point
/// @return The elapsed time between now and *begin* in seconds
auto elapsed(Time const& begin) -> double;

/// Used for measuring elapsed time since object creation.
class Stopwatch
{
public:
    /// Construct a Stopwatch object and start measuring time.
    Stopwatch();

    /// Start measuring time.
    auto start() -> void;

    /// Pause measuring time.
    auto pause() -> void;

    /// Reset the stopwatch.
    auto reset() -> void;

    /// Get the accumulated elapsed time (in seconds) between calls to methods @ref start and @ref pause.
    auto time() const -> double;

private:
    /// The auxiliary time variable marking the start of timing.
    Time mstart;

    /// The elapsed time measured by the stopwatch.
    double melapsed = 0.0;
};

} // namespace Reaktoro
