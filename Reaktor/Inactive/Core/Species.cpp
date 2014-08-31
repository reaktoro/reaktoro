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

#include "Species.hpp"


// C++ includes
#include <sstream>

// Reaktor includes
#include <Reaktor/Utils/ElementUtils.hpp>
#include <Reaktor/Utils/StringUtils.hpp>

namespace Reaktor {
namespace internal {

auto defaultChemicalPotentialFn() -> ChemicalPotentialFn
{
    return [](double, double) { return 0.0; };
}

auto defaultDensity() -> DensityFn
{
    return [](double, double) { return 0.0; };
}

} /* namespace internal */

using namespace internal;

Species::Species()
: charge$(0), type$("none"),
  chemical_potential$(defaultChemicalPotentialFn()),
  density$(defaultDensity())
{}

auto Species::setName(const std::string& name) -> void
{
	name$ = name;
}

auto Species::setFormula(const std::string& formula) -> void
{
    formula$ = formula;
}

auto Species::setElements(const std::map<std::string, double>& elements) -> void
{
    elements$ = elements;

    setMolarMass(Reaktor::molarMass(elements));
}

auto Species::setElements(const std::string& elements) -> void
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

auto Species::setMolarMass(units::MolarMass value) -> void
{
	molar_mass$ = value;
}

auto Species::setCharge(double value) -> void
{
	charge$ = value;
}

auto Species::setType(const std::string& type) -> void
{
    type$ = type;
}

auto Species::setChemicalPotential(const ChemicalPotentialFn& function) -> void
{
	chemical_potential$ = function;
}

auto Species::setChemicalPotential(double constant) -> void
{
    chemical_potential$ = [=](double,double){ return constant; };
}

auto Species::setDensity(const DensityFn& density) -> void
{
    density$ = density;
}

auto Species::setDensity(double density) -> void
{
    auto density_func = [=](double, double) { return density; };

    setDensity(density_func);
}

auto Species::name() const -> const std::string&
{
	return name$;
}

auto Species::formula() const -> std::string
{
    return formula$;
}

auto Species::elements() const -> const std::map<std::string, double>&
{
	return elements$;
}

auto Species::molarMass() const -> units::MolarMass
{
	return molar_mass$;
}

auto Species::charge() const -> double
{
	return charge$;
}

auto Species::type() const -> const std::string&
{
    return type$;
}

auto Species::containsElement(const std::string& element) const -> bool
{
    return elements$.count(element);
}

auto Species::elementAtoms(const std::string& element) const -> double
{
    auto iter = elements$.find(element);
    return iter != elements$.end() ? iter->second : 0.0;
}

auto Species::elementAtoms() const -> std::vector<double>
{
    std::vector<double> atoms;
    atoms.reserve(elements$.size());
    for(const auto& entry : elements$)
        atoms.push_back(entry.second);
    return atoms;
}

auto Species::elementNames() const -> std::vector<std::string>
{
    std::vector<std::string> names;
    names.reserve(elements$.size());
    for(const auto& entry : elements$)
        names.push_back(entry.first);
    return names;
}

auto Species::chemicalPotential() const -> const ChemicalPotentialFn&
{
	return chemical_potential$;
}

auto Species::chemicalPotential(double T, double P) const -> double
{
	return chemical_potential$(T, P);
}

auto Species::density(double T, double P) const -> double
{
    return density$(T, P);
}

auto Species::molarDensity(double T, double P) const -> double
{
    return density(T, P)/molarMass().in(unit(kg)/unit(mol));
}

auto Species::molarVolume(double T, double P) const -> double
{
    return 1.0/molarDensity(T, P);
}

auto Species::specificVolume(double T, double P) const -> double
{
    return 1.0/density(T, P);
}

Species::operator std::string() const
{
	return name$;
}

auto names(const std::vector<Species>& species) -> std::vector<std::string>
{
    std::vector<std::string> names(species.size());
    for(unsigned i = 0; i < species.size(); ++i)
        names[i] = species[i].name();
    return names;
}

auto charges(const std::vector<Species>& species) -> Vector
{
    Vector charges(species.size());
    for(unsigned i = 0; i < species.size(); ++i)
        charges[i] = species[i].charge();
    return charges;
}

auto operator<<(std::ostream& out, const Species& species) -> std::ostream&
{
    auto str = [](const std::map<std::string, double>& elements)
    {
        std::stringstream ss;
        for(const auto& pair : elements)
            ss << pair.first << "(" << pair.second << ")";
        return ss.str();
    };

    out << "Species: " << species.name() << std::endl;
    out << " Chemical Formula: " << species.formula() << std::endl;
    out << " Elemental Formula: " << str(species.elements()) << std::endl;
    out << " MolarMass: " << species.molarMass() << std::endl;

    return out;
}

} // namespace Reaktor
