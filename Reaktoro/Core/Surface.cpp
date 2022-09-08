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
    String name = "NotSpecified";

    /// The names of the phases composing this surface.
    Pair<String, String> phases = {"NotSpecified", "NotSpecified"};

    /// The indices of the phases composing this surface.
    Pair<Index, Index> iphases = {-1, -1};

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

Surface::Surface(String const& name, String const& phase1, Index iphase1, String const& phase2, Index iphase2)
: Surface()
{
    pimpl->name = name;
    pimpl->phases = {phase1, phase2};
    pimpl->iphases = {iphase1, iphase2};
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

auto Surface::withPhaseNames(String const& phase1, String const& phase2) const -> Surface
{
    Surface copy = clone();
    copy.pimpl->phases = {phase1, phase2};
    return copy;
}

auto Surface::withPhaseIndices(Index iphase1, Index iphase2) const -> Surface
{
    Surface copy = clone();
    copy.pimpl->iphases = {iphase1, iphase2};
    return copy;
}

auto Surface::name() const -> String
{
    return pimpl->name;
}

auto Surface::phaseNames() const -> Pair<String, String> const&
{
    return pimpl->phases;
}

auto Surface::phaseIndices() const -> Pair<Index, Index> const&
{
    return pimpl->iphases;
}

auto Surface::equivalent(Surface const& another) const -> bool
{
    auto const& [phase1, phase2] = phaseNames();
    auto const& [other1, other2] = another.phaseNames();

    auto const& [iphase1, iphase2] = phaseIndices();
    auto const& [iother1, iother2] = another.phaseIndices();

    auto const same_phase_names = (phase1 == other1 && phase2 == other2) || (phase1 == other2 && phase2 == other1);
    auto const same_phase_indcs = (iphase1 == iother1 && iphase2 == iother2) || (iphase1 == iother2 && iphase2 == iother1);

    return same_phase_names && same_phase_indcs;
}

auto Surface::equivalent(String const& otherphase1, String const& otherphase2) const -> bool
{
    auto const& [phase1, phase2] = phaseNames();
    return (phase1 == otherphase1 && phase2 == otherphase2) || (phase1 == otherphase2 && phase2 == otherphase1);
}

auto Surface::equivalent(Index iotherphase1, Index iotherphase2) const -> bool
{
    auto const& [iphase1, iphase2] = phaseIndices();
    return (iphase1 == iotherphase1 && iphase2 == iotherphase2) || (iphase1 == iotherphase2 && iphase2 == iotherphase1);
}

auto Surface::equivalent(StringOrIndex const& phase1, StringOrIndex const& phase2) const -> bool
{
    const auto [nphase1, nphase2] = phaseNames();
    const auto [iphase1, iphase2] = phaseIndices();

    auto same = [](StringOrIndex const& si, String const& s, Index i)
    {
        errorif(si.valueless_by_exception(), "Invalid StringOrIndex object.");
        if(auto pval = std::get_if<Index>(&si))  return *pval == i;
        if(auto pval = std::get_if<int>(&si))    return *pval == i;
        if(auto pval = std::get_if<String>(&si)) return *pval == s;
        return false;
    };

    return
        same(phase1, nphase1, iphase1) && same(phase2, nphase2, iphase2) ||
        same(phase2, nphase1, iphase1) && same(phase1, nphase2, iphase2);
}

auto operator<(const Surface& lhs, const Surface& rhs) -> bool
{
    return lhs.name() < rhs.name();
}

auto operator==(const Surface& lhs, const Surface& rhs) -> bool
{
    return lhs.name() == rhs.name() && lhs.phaseNames() == rhs.phaseNames() && lhs.phaseIndices() == rhs.phaseIndices();
}

} // namespace Reaktoro
