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

// Boost includes
#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

using namespace boost::python;

// Reaktor includes
#include <Reaktor/Core/Species.hpp>

// PyReator includes
#include <python/Utils/Converters.hpp>

namespace Reaktor {

auto createSpecies(std::string name, std::vector<Element> elements, std::vector<double> atoms, double charge) -> boost::shared_ptr<Species>
{
	SpeciesData data;
	data.name = std::move(name);
	data.elements = std::move(elements);
	data.atoms = std::move(atoms);
	data.charge = std::move(charge);
	return boost::make_shared<Species>(data);
}

auto exportSpecies() -> void
{
	class_<SpeciesData>("SpeciesData")
		.def_readwrite("name", &SpeciesData::name)
		.def_readwrite("formula", &SpeciesData::formula)
		.def_readwrite("elements", &SpeciesData::elements)
		.def_readwrite("atoms", &SpeciesData::atoms)
		.def_readwrite("charge", &SpeciesData::charge)
		.def_readwrite("molar_mass", &SpeciesData::molar_mass)
		;

	class_<Species>("Species")
		.def(init<>())
		.def(init<const SpeciesData&>())
		.def("__init__", make_constructor(createSpecies, default_call_policies(),
			(arg("name"), arg("elements"), arg("atoms"), arg("charge"))))
		.def("name", &Species::name, return_value_policy<copy_const_reference>())
		.def("formula", &Species::formula, return_value_policy<copy_const_reference>())
		.def("elements", &Species::elements, return_value_policy<copy_const_reference>())
		.def("atoms", &Species::atoms, return_value_policy<copy_const_reference>())
		.def("charge", &Species::charge)
		.def("molarMass", &Species::molarMass)
		;

	export_std_vector<Species>("SpeciesVector");

	def("atoms", atoms);
	def("collectElements", collectElements);
	def("collectCharges", collectCharges);
	def("collectMolarMasses", collectMolarMasses);
}

} // namespace Reaktor
