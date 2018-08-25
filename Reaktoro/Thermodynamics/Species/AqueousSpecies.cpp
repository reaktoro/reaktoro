// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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


