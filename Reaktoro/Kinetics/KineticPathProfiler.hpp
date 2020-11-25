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
#include <deque>
#include <ostream>
#include <memory>

// Reaktoro includes
#include <Reaktoro/Kinetics/KineticPathAnalysis.hpp>

namespace Reaktoro {

/// Provide mechanisms for analysing the operations in a kinetic path simulation.
class KineticPathProfiler
{
public:
    /// Construct a default instance of KineticPathProfiler.
    KineticPathProfiler();

    /// Construct a copy of a KineticPathProfiler instance.
    KineticPathProfiler(const KineticPathProfiler& other);

    /// Destroy this KineticPathProfiler instance.
    virtual ~KineticPathProfiler();

    /// Assign a copy of an KineticPathProfiler instance.
    auto operator=(KineticPathProfiler other) -> KineticPathProfiler&;

    /// Update the profiler with the result of the last kinetic path time step.
    auto update(const KineticResult& result) -> void;

    /// Return the performance analysis of all operations in the kinetic path simulation.
    auto analysis() const -> KineticPathAnalysis;

    /// Return all collected results of the kinetic path simulation.
    auto results() const -> const std::deque<KineticResult>&;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

auto operator<<(std::ostream& out, const KineticPathProfiler& prof) -> std::ostream&;

} // namespace Reaktoro
