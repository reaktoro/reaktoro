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

namespace Reaktor {

struct Species::Impl
{
	Impl()
	: molar_mass(0.0), charge(0.0)
	{}

	/// The name of the chemical species
    std::string name;

    /// The formula of the chemical species
    std::string formula;

    /// The elements that compose the chemical species
    std::vector<std::string> elements;

    /// The coefficients of the elements that compose the chemical species
    std::vector<double> coefficients;

    /// The molar mass of the chemical species
    double molar_mass;

    /// The electrical charge of the chemical species
    double charge;

    /// The standard chemical potential function of the chemical species
    ChemicalPotential chemical_potential;
};

Species::Species()
: pimpl(new Impl())
{}

Species::Species(const Species& other)
: pimpl(new Impl(*other.pimpl))
{}

Species::~Species()
{}

auto Species::operator=(Species other) -> Species&
{
	pimpl = std::move(other.pimpl);
	return *this;
}

auto Species::setName(const std::string& name) -> Species&
{
	pimpl->name = name;
    return *this;
}

auto Species::setFormula(const std::string& formula) -> Species&
{
    pimpl->formula = formula;
    return *this;
}

auto Species::setElements(const std::vector<std::string>& elements) -> Species&
{
	pimpl->elements = elements;
    return *this;
}

auto Species::setCoefficients(const std::vector<double>& coefficients) -> Species&
{
	pimpl->coefficients = coefficients;
    return *this;
}

auto Species::setMolarMass(double value) -> Species&
{
	pimpl->molar_mass = value;
    return *this;
}

auto Species::setCharge(double value) -> Species&
{
	pimpl->charge = value;
    return *this;
}

auto Species::setChemicalPotential(const ChemicalPotential& function) -> Species&
{
	pimpl->chemical_potential = function;
    return *this;
}

auto Species::name() const -> const std::string&
{
	return pimpl->name;
}

auto Species::formula() const -> const std::string&
{
	return pimpl->formula;
}

auto Species::elements() const -> const std::vector<std::string>&
{
	return pimpl->elements;
}

auto Species::coefficients() const -> const std::vector<double>&
{
	return pimpl->coefficients;
}

auto Species::molarMass() const -> double
{
	return pimpl->molar_mass;
}

auto Species::charge() const -> double
{
	return pimpl->charge;
}

auto Species::chemicalPotential() const -> const ChemicalPotential&
{
	return pimpl->chemical_potential;
}

auto operator<<(std::ostream& out, const Species& species) -> std::ostream&
{
	// todo improve operator<< of Species
    out << "Species" << std::endl;
    out << "    name:         " << species.name() << std::endl;
    out << "    formula:      " << species.formula() << std::endl;
//    out << "    elements:     " << species.elements() << std::endl;
//    out << "    coefficients: " << species.coefficients() << std::endl;
    out << "    charge:       " << species.charge() << std::endl;
    out << "    molar-mass:   " << species.molarMass() << " mol/kg" << std::endl;

    return out;
}

} /* namespace Reaktor */
