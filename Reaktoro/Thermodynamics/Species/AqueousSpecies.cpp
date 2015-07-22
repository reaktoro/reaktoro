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

#include "AqueousSpecies.hpp"

namespace Reaktoro {

struct AqueousSpecies::Impl
{
    /// The electrical charge of the aqueous species
    double charge;

    /// The dissociation of a neutral aqueous species into charged species.
    std::map<std::string, double> dissociation;

    /// The thermodynamic data of the aqueous species.
    AqueousSpeciesThermoData thermo;
};

AqueousSpecies::AqueousSpecies()
: pimpl(new Impl())
{}

AqueousSpecies::AqueousSpecies(const Species& species)
: Species(species), pimpl(new Impl())
{}

AqueousSpecies::AqueousSpecies(const AqueousSpecies& other)
: Species(other), pimpl(new Impl(*other.pimpl))
{}

AqueousSpecies::~AqueousSpecies()
{}

auto AqueousSpecies::operator=(AqueousSpecies other) -> AqueousSpecies&
{
    Species::operator=(other);
    pimpl = std::move(other.pimpl);
    return *this;
}

auto AqueousSpecies::setCharge(double value) -> void
{
    pimpl->charge = value;
}

auto AqueousSpecies::setDissociation(const std::map<std::string, double>& dissociation) -> void
{
    pimpl->dissociation = dissociation;
}

auto AqueousSpecies::setThermoData(const AqueousSpeciesThermoData& thermo) -> void
{
    pimpl->thermo = thermo;
}

auto AqueousSpecies::charge() const -> double
{
    return pimpl->charge;
}

auto AqueousSpecies::dissociation() const -> const std::map<std::string, double>&
{
    return pimpl->dissociation;
}

auto AqueousSpecies::thermoData() const -> const AqueousSpeciesThermoData&
{
    return pimpl->thermo;
}

} // namespace Reaktoro


