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

// Reaktor includes
#include <Reaktor/Core/Utils.hpp>

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

auto Species::elements() const -> const std::vector<Element>&
{
    return data->elements;
}

auto Species::atoms() const -> const std::vector<double>&
{
    return data->atoms;
}

auto Species::charge() const -> double
{
	return data->charge;
}

auto Species::molarMass() const -> double
{
    return data->molar_mass;
}

auto operator<(const Species& lhs, const Species& rhs) -> bool
{
	return lhs.name() < rhs.name();
}

auto operator==(const Species& lhs, const Species& rhs) -> bool
{
	return lhs.name() == rhs.name();
}

auto atoms(const Element& element, const Species& species) -> double
{
    Index idx = index(element, species.elements());
    return idx < species.elements().size() ? species.atoms()[idx] : 0.0;
}

auto collectElements(const std::vector<Species>& species) -> std::vector<Element>
{
    std::set<Element> elements;
    for(const Species& iter : species)
        elements.insert(iter.elements().begin(), iter.elements().end());
    return std::vector<Element>(elements.begin(), elements.end());
}

auto collectCharges(const std::vector<Species>& species) -> Vector
{
    Vector charges(species.size());
    for(unsigned i = 0; i < species.size(); ++i)
        charges[i] = species[i].charge();
    return charges;
}

auto collectMolarMasses(const std::vector<Species>& species) -> Vector
{
    Vector molar_masses(species.size());
    for(unsigned i = 0; i < species.size(); ++i)
        molar_masses[i] = species[i].molarMass();
    return molar_masses;
}

} // namespace Reaktor
