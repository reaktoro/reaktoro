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

#include "Surfaces.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>

namespace Reaktoro {

Surfaces::Surfaces()
{
}

auto Surfaces::add(String const& phase1, String const& phase2) -> void
{
    errorif(phase1.empty() || phase2.empty(), "Expecting a non-empty phase name when registering a surface.");
    m_surfaces.push_back(Surface(phase1 + ":" + phase2).withPhases(phase1, phase2));
}

auto Surfaces::add(String const& phase) -> void
{
    errorif(phase.empty(), "Expecting a non-empty phase name when registering a surface.");
    m_surfaces.push_back(Surface(phase).withPhases(phase, phase));
}

auto Surfaces::surfaces() const -> Vec<Surface> const&
{
    return m_surfaces;
}

} // namespace Reaktoro
