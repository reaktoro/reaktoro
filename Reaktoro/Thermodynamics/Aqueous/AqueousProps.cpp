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

#include "AqueousProps.hpp"

// C++ includes

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>

namespace Reaktoro {

struct AqueousProps::Impl
{
    /// The underlying Phase object for the aqueous phase.
    Phase phase;

    /// Construct a default AqueousProps::Impl object.
    Impl(const Phase& phase)
    : phase(phase)
    {}
};

AqueousProps::AqueousProps(const Phase& phase)
: PhaseChemicalProps(phase), pimpl(new Impl(phase))
{
    error(true, "AqueousProps class has not been implemented yet.");
}

AqueousProps::AqueousProps(const AqueousProps& other)
: PhaseChemicalProps(other.phase()), pimpl(new Impl(*other.pimpl))
{
    error(true, "AqueousProps class has not been implemented yet.");
}

AqueousProps::~AqueousProps()
{}

auto AqueousProps::operator=(AqueousProps other) -> AqueousProps&
{
    PhaseChemicalProps::operator=(other);
    pimpl = std::move(other.pimpl);
    return *this;
}

auto AqueousProps::ionicStrength() const -> real
{
    error(true, "AqueousProps class has not been implemented yet.");
}

auto AqueousProps::pH() const -> real
{
    error(true, "AqueousProps class has not been implemented yet.");
}

auto AqueousProps::pE() const -> real
{
    error(true, "AqueousProps class has not been implemented yet.");
}

auto AqueousProps::Eh() const -> real
{
    error(true, "AqueousProps class has not been implemented yet.");
}

auto AqueousProps::alkalinity() const -> real
{
    error(true, "AqueousProps class has not been implemented yet.");
}

auto AqueousProps::phase() const -> const Phase&
{
    return pimpl->phase;
}

} // namespace Reaktoro
