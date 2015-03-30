// Reaktor is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "GaseousSpecies.hpp"

namespace Reaktor {

struct GaseousSpecies::Impl
{
    /// The thermodynamic data of the gaseous species.
    GaseousSpeciesThermoData thermo;
};

GaseousSpecies::GaseousSpecies()
: pimpl(new Impl())
{}

GaseousSpecies::GaseousSpecies(const GeneralSpecies& species)
: GeneralSpecies(species), pimpl(new Impl())
{}

GaseousSpecies::GaseousSpecies(const GaseousSpecies& other)
: GeneralSpecies(other), pimpl(new Impl(*other.pimpl))
{}

GaseousSpecies::~GaseousSpecies()
{}

auto GaseousSpecies::operator=(GaseousSpecies other) -> GaseousSpecies&
{
    GeneralSpecies::operator=(other);
    pimpl = std::move(other.pimpl);
    return *this;
}

auto GaseousSpecies::setThermoData(const GaseousSpeciesThermoData& thermo) -> void
{
    pimpl->thermo = thermo;
}

auto GaseousSpecies::thermoData() const -> const GaseousSpeciesThermoData&
{
    return pimpl->thermo;
}

} // namespace Reaktor
