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
#include <Reaktoro/Core/Surface.hpp>

namespace Reaktoro {

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

    /// Return the Surface objects already constructed.
    auto surfaces() const -> Vec<Surface> const&;

private:
    /// The registered surfaces.
    Vec<Surface> m_surfaces;
};

} // namespace Reaktoro
