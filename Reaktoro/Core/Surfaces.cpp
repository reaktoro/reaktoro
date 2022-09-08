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
#include <Reaktoro/Core/PhaseList.hpp>
#include <Reaktoro/Core/Surface.hpp>

namespace Reaktoro {

Surfaces::Surfaces()
{
}

auto Surfaces::add(String const& phase1, String const& phase2) -> void
{
    errorif(phase1.empty() || phase2.empty(), "Expecting a non-empty phase name when registering a surface.");
    surfaces.push_back({ phase1, phase2});
}

auto Surfaces::add(String const& phase) -> void
{
    errorif(phase.empty(), "Expecting a non-empty phase name when registering a surface.");
    surfaces.push_back({ phase, phase});
}

auto Surfaces::data() const -> Pairs<String, String> const&
{
    return surfaces;
}

auto Surfaces::convert(PhaseList const& phases) const -> Vec<Surface>
{
    auto createSurface = [&](auto phasepair)
    {
        auto const nphase1 = phasepair.first;
        auto const nphase2 = phasepair.second;
        auto const iphase1 = phases.findWithName(nphase1);
        auto const iphase2 = phases.findWithName(nphase2);
        errorif(iphase1 >= phases.size(), "Expecting a name of a phase that exist in the list of assembled phase, but got instead `", nphase1, "`.");
        errorif(iphase2 >= phases.size(), "Expecting a name of a phase that exist in the list of assembled phase, but got instead `", nphase2, "`.");
        auto const id = nphase1 != nphase2 ? nphase1 + ":" + nphase2 : nphase1;
        return Surface(id, nphase1, iphase1, nphase2, iphase2);
    };

    return vectorize(surfaces, RKT_LAMBDA(x, createSurface(x)));
}

} // namespace Reaktoro
