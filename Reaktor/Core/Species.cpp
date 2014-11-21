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

#include "Species.hpp"

// C++ includes
#include <set>

namespace Reaktor {

Species::Species()
: data(new SpeciesData())
{}

Species::Species(const SpeciesData& data)
: data(new SpeciesData(data))
{}

auto Species::name() const -> const std::string&
{
	return data->name;
}

auto Species::formula() const -> const std::string&
{
	return data->formula;
}

auto Species::components() const -> const ComponentList&
{
    return data->components;
}

auto Species::stoichiometries() const -> const std::vector<double>&
{
    return data->stoichiometries;
}

auto Species::charge() const -> double
{
	return data->charge;
}

auto Species::molarMass() const -> double
{
    return data->molar_mass;
}

auto components(const SpeciesList& species) -> ComponentList
{
    std::set<Component> elements;
    for(const Species& iter : species)
        elements.insert(iter.components().begin(), iter.components().end());
    return ComponentList(elements.begin(), elements.end());
}

auto charges(const SpeciesList& species) -> std::vector<double>
{
    std::vector<double> charges;
    charges.reserve(species.size());
    for(const Species& iter : species)
        charges.push_back(iter.charge());
    return charges;
}

auto molarMasses(const SpeciesList& species) -> std::vector<double>
{
    std::vector<double> molar_masses;
    molar_masses.reserve(species.size());
    for(const Species& iter : species)
        molar_masses.push_back(iter.molarMass());
    return molar_masses;
}

} // namespace Reaktor
