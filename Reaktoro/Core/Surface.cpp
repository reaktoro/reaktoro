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

    /// The names of the phases composing this surface.
    Pair<String, String> phases;

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

Surface::Surface(String const& name, String const& phase1, String const& phase2)
: Surface()
{
    pimpl->name = name;
    pimpl->phases = {phase1, phase2};
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

auto Surface::withPhases(String const& phase1, String const& phase2) const -> Surface
{
    Surface copy = clone();
    copy.pimpl->phases = {phase1, phase2};
    return copy;
}

auto Surface::name() const -> String
{
    return pimpl->name;
}

auto Surface::phases() const -> Pair<String, String> const&
{
    return pimpl->phases;
}

auto Surface::equivalent(Surface const& another) const -> bool
{
    auto const& [phase1, phase2] = phases();
    auto const& [other1, other2] = another.phases();
    return (phase1 == other1 && phase2 == other2) || (phase1 == other2 && phase2 == other1);
}

auto Surface::equivalent(String const& otherphase1, String const& otherphase2) const -> bool
{
    auto const& [phase1, phase2] = phases();
    return (phase1 == otherphase1 && phase2 == otherphase2) || (phase1 == otherphase2 && phase2 == otherphase1);
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
