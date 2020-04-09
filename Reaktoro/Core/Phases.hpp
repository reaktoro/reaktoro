// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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

// Reaktoro includes
#include <Reaktoro/Core/PhaseEditor.hpp>

namespace Reaktoro {

/// The base type for all other classes defining phase types.
/// @ingroup Core
class Phases
{
public:
    /// Construct a Phases object.
    template<typename... Args>
    explicit Phases(const ThermoEngine& engine, Args&&... phases);

private:
    std::vector<std::unique_ptr<PhaseEditor>> phases;
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

/// Compare two Phases instances for less than
auto operator<(const Phases& lhs, const Phases& rhs) -> bool;

/// Compare two Phases instances for equality
auto operator==(const Phases& lhs, const Phases& rhs) -> bool;

} // namespace Reaktoro
