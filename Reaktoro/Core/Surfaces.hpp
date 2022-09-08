// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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
#include <Reaktoro/Common/Types.hpp>

namespace Reaktoro {

// Forward declarations
class PhaseList;
class Surface;

/// Used to represent the interface surface between two phases in a chemical system.
class Surfaces
{
public:
    /// Construct a Surfaces object.
    Surfaces();

    /// Register a new surface between phases with names `phase1` and `phase2`.
    auto add(String const& phase1, String const& phase2) -> void;

    /// Register a new surface surrounding phase with name `phase`.
    auto add(String const& phase) -> void;

    /// Return the pairs of phase names already registered.
    auto data() const -> Pairs<String, String> const&;

    /// Convert this Surfaces object into a vector of Surface objects.
    auto convert(PhaseList const& phases) const -> Vec<Surface>;

private:
    /// The registered surfaces.
    Pairs<String, String> surfaces;
};

} // namespace Reaktoro
