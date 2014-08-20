/*
 * Reaktor is a C++ library for computational reaction modelling.
 *
 * Copyright (C) 2014 Allan Leal
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "GeneralSpecies.hpp"


// C++ includes
#include <sstream>

// Reaktor includes
#include <Reaktor/Utils/ElementUtils.hpp>
#include <Reaktor/Utils/StringUtils.hpp>

namespace Reaktor {

GeneralSpecies::GeneralSpecies()
: molar_mass$(0), charge$(0)
{}

auto GeneralSpecies::setName(const std::string& name) -> void
{
	name$ = name;
}

auto GeneralSpecies::setFormula(const std::string& formula) -> void
{
    formula$ = formula;
}

auto GeneralSpecies::setElements(const std::map<std::string, double>& elements) -> void
{
	elements$ = elements;

	setMolarMass(Reaktor::molarMass(elements));
}

auto GeneralSpecies::setElements(const std::string& elements) -> void
{
    // The map with the element names and number of atoms
    std::map<std::string, double> element_map;

    // Split the elemental formula in words delimited by ( or )
    auto words = split(elements, "()");

    // Create the pair entries
    for(unsigned i = 0; i < words.size(); i += 2)
        element_map.insert({words[i], tofloat(words[i + 1])});

    // Set the elements of the species
    setElements(element_map);
}

auto GeneralSpecies::setMolarMass(units::MolarMass value) -> void
{
	molar_mass$ = value;
}

auto GeneralSpecies::setCharge(double value) -> void
{
	charge$ = value;
}

auto GeneralSpecies::name() const -> const std::string&
{
	return name$;
}

auto GeneralSpecies::formula() const -> std::string
{
    return formula$;
}

auto GeneralSpecies::elements() const -> const std::map<std::string, double>&
{
	return elements$;
}

auto GeneralSpecies::molarMass() const -> units::MolarMass
{
	return molar_mass$;
}

auto GeneralSpecies::charge() const -> double
{
	return charge$;
}

auto GeneralSpecies::containsElement(const std::string& element) const -> bool
{
    return elements$.count(element);
}

auto GeneralSpecies::elementAtoms(const std::string& element) const -> double
{
    auto iter = elements$.find(element);
    return iter != elements$.end() ? iter->second : 0.0;
}

auto GeneralSpecies::elementAtoms() const -> std::vector<double>
{
    std::vector<double> atoms;
    atoms.reserve(elements$.size());
    for(const auto& entry : elements$)
        atoms.push_back(entry.second);
    return atoms;
}

auto GeneralSpecies::elementNames() const -> std::vector<std::string>
{
    std::vector<std::string> names;
    names.reserve(elements$.size());
    for(const auto& entry : elements$)
        names.push_back(entry.first);
    return names;
}

GeneralSpecies::operator std::string() const
{
	return name$;
}

auto operator<<(std::ostream& out, const GeneralSpecies& species) -> std::ostream&
{
    auto str = [](const std::map<std::string, double>& elements)
    {
        std::stringstream ss;
        for(const auto& pair : elements)
            ss << pair.first << "(" << pair.second << ")";
        return ss.str();
    };

    out << "Species: " << species.name() << std::endl;
    out << "  chemical formula: " << species.formula() << std::endl;
    out << "  elemental formula: " << str(species.elements()) << std::endl;
    out << "  molar mass: " << species.molarMass() << " g/mol" << std::endl;

    return out;
}

} /* namespace Reaktor */
