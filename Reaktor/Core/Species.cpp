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
    /// The name of the species
	std::string name;

	/// The chemical formula of the species
	std::string formula;

	/// The names of the elements that compose the species
    std::vector<std::string> elements;

    /// The number of atoms of the elements that compose the species
    std::vector<double> element_atoms;

	/// The molar mass of the species (in units of kg/mol)
	double molar_mass;

	/// The electrical charge of the species
	double charge;

	/// The thermodynamic model of the species
	SpeciesThermoModel thermo_model;
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

auto Species::setElementAtoms(const std::vector<double>& element_atoms) -> Species&
{
    pimpl->element_atoms = element_atoms;
    return *this;
}

auto Species::setMolarMass(double val) -> Species&
{
	pimpl->molar_mass = val;
	return *this;
}

auto Species::setCharge(double val) -> Species&
{
	pimpl->charge = val;
	return *this;
}

auto Species::setThermoModel(const SpeciesThermoModel& thermo_model) -> Species&
{
	pimpl->thermo_model = thermo_model;
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

auto Species::elementAtoms() const -> const std::vector<double>&
{
    return pimpl->element_atoms;
}

auto Species::molarMass() const -> double
{
	return pimpl->molar_mass;
}

auto Species::charge() const -> double
{
	return pimpl->charge;
}

auto Species::thermoModel() const -> const SpeciesThermoModel&
{
	return pimpl->thermo_model;
}

} // namespace Reaktor
