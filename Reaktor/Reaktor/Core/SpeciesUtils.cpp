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

#include "SpeciesUtils.hpp"

// C++ includes
#include <algorithm>

// Reaktor includes
#include <Reaktor/Core/Species.hpp>

namespace Reaktor {

auto numElements(const Species& species) -> unsigned
{
	return species.elements().size();
}

auto indexElement(const Species& species, const std::string& element) -> Index
{
	const auto& begin = species.elements().begin();
	const auto& end = species.elements().end();
	return std::find(begin, end, element) - begin;
}

auto containsElement(const Species& species, const std::string& element) -> bool
{
    return indexElement(species, element) < numElements(species);
}

auto elementAtoms(const Species& species, const std::string& element) -> double
{
	const Index ielement = indexElement(species, element);
    return ielement < numElements(species) ? species.coefficients()[ielement] : 0.0;
}

auto chemicalPotential(const Species& species, double T, double P) -> double
{
	return species.chemicalPotential()(T, P);
}

auto names(const std::vector<Species>& species) -> std::vector<std::string>
{
    std::vector<std::string> names(species.size());
    for(unsigned i = 0; i < species.size(); ++i)
        names[i] = species[i].name();
    return names;
}

auto charges(const std::vector<Species>& species) -> std::vector<double>
{
    std::vector<double> charges(species.size());
    for(unsigned i = 0; i < species.size(); ++i)
        charges[i] = species[i].charge();
    return charges;
}

} /* namespace Reaktor */
