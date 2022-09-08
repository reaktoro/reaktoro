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

#include "Surface.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>

namespace Reaktoro {

struct Surface::Impl
{
    /// The unique name of the surface.
    String name;

    /// The 1st phase composing this surface.
    Phase phase1;

    /// The 2nd phase composing this surface.
    Phase phase2;

    /// Construct a default Surface::Impl object.
    Impl()
    {}
};

Surface::Surface()
: pimpl(new Impl())
{}

Surface::Surface(String const& name)
: Surface()
{
    pimpl->name = name;
}

auto Surface::clone() const -> Surface
{
    Surface surface;
    *surface.pimpl = *pimpl;
    return surface;
}

auto Surface::withName(String const& name) const -> Surface
{
    Surface copy = clone();
    copy.pimpl->name = name;
    return copy;
}

auto Surface::withPhases(Phase const& phase1, Phase const& phase2) const -> Surface
{
    Surface copy = clone();
    copy.pimpl->phase1 = phase1;
    copy.pimpl->phase2 = phase2;
    return copy;
}

auto Surface::name() const -> String
{
    return pimpl->name;
}

auto Surface::phases() const -> Pair<Phase const&, Phase const&>
{
    return {pimpl->phase1, pimpl->phase2};
}

auto operator<(const Surface& lhs, const Surface& rhs) -> bool
{
    return lhs.name() < rhs.name();
}

auto operator==(const Surface& lhs, const Surface& rhs) -> bool
{
    return lhs.name() == rhs.name() && lhs.phases() == rhs.phases();
}

} // namespace Reaktoro
