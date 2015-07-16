// Reaktoro is a C++ library for computational reaction modelling.
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

#include "MineralSpecies.hpp"

namespace Reaktoro {

struct MineralSpecies::Impl
{
    /// The thermodynamic data of the mineral species.
    MineralSpeciesThermoData thermo;
};

MineralSpecies::MineralSpecies()
: pimpl(new Impl())
{}

MineralSpecies::MineralSpecies(const Species& species)
: Species(species), pimpl(new Impl())
{}

MineralSpecies::MineralSpecies(const MineralSpecies& other)
: Species(other), pimpl(new Impl(*other.pimpl))
{}

MineralSpecies::~MineralSpecies()
{}

auto MineralSpecies::operator=(MineralSpecies other) -> MineralSpecies&
{
    Species::operator=(other);
    pimpl = std::move(other.pimpl);
    return *this;
}

auto MineralSpecies::setThermoData(const MineralSpeciesThermoData& thermo) -> void
{
    pimpl->thermo = thermo;
}

auto MineralSpecies::thermoData() const -> const MineralSpeciesThermoData&
{
    return pimpl->thermo;
}

} // namespace Reaktoro
