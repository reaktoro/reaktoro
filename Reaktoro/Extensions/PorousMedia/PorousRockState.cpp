// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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

#include "PorousRockState.hpp"

// C++ includes

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>

namespace Reaktoro {

struct PorousRockState::Impl
{
    Impl(const ChemicalSystem& system)
    {
    }
};

PorousRockState::PorousRockState(const ChemicalSystem& system)
: ChemicalState(system), pimpl(new Impl(system))
{}

PorousRockState::PorousRockState(const PorousRockState& other)
: ChemicalState(other.system()), pimpl(new Impl(*other.pimpl))
{}

PorousRockState::~PorousRockState()
{}

auto PorousRockState::operator=(PorousRockState other) -> PorousRockState&
{
    ChemicalState::operator=(other);
    pimpl = std::move(other.pimpl);
    return *this;
}

auto PorousRockState::setPorosity(real value) -> void
{
    error(true, "PorousRockState has not been implemented yet.");
}

auto PorousRockState::setVolumeFractionAmongFluids(Index iphase, real value) -> void
{
    error(true, "PorousRockState has not been implemented yet.");
}

auto PorousRockState::setVolumeFractionAmongFluids(String phasename, real value) -> void
{
    error(true, "PorousRockState has not been implemented yet.");
}

auto PorousRockState::setVolumeFractionAmongSolids(Index iphase, real value) -> void
{
    error(true, "PorousRockState has not been implemented yet.");
}

auto PorousRockState::setVolumeFractionAmongSolids(String phasename, real value) -> void
{
    error(true, "PorousRockState has not been implemented yet.");
}

auto PorousRockState::equilibrate() -> void
{
    error(true, "PorousRockState has not been implemented yet.");
}

} // namespace Reaktoro
