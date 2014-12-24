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

#include "PySpecies.hpp"

// Boost includes
#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
namespace py = boost::python;

// Reaktor includes
#include <Reaktor/Reaktor.hpp>

// PyReator includes
#include <PyReaktor/Utils/Converters.hpp>

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

auto export_Species() -> void
{
	py::class_<SpeciesData>("SpeciesData")
		.def_readwrite("name", &SpeciesData::name)
		.def_readwrite("formula", &SpeciesData::formula)
		.def_readwrite("elements", &SpeciesData::elements)
		.def_readwrite("atoms", &SpeciesData::atoms)
		.def_readwrite("charge", &SpeciesData::charge)
		.def_readwrite("molar_mass", &SpeciesData::molar_mass)
		;

	py::class_<Species>("Species")
		.def(py::init<>())
		.def(py::init<const SpeciesData&>())
		.def("__init__", py::make_constructor(createSpecies, py::default_call_policies(),
			(py::arg("name"), py::arg("elements"), py::arg("atoms"), py::arg("charge"))))
		.def("name", &Species::name, py::return_value_policy<py::copy_const_reference>())
		.def("formula", &Species::formula, py::return_value_policy<py::copy_const_reference>())
		.def("elements", &Species::elements, py::return_value_policy<py::copy_const_reference>())
		.def("atoms", &Species::atoms, py::return_value_policy<py::copy_const_reference>())
		.def("charge", &Species::charge)
		.def("molarMass", &Species::molarMass)
		;

	py::def("atoms", atoms);
	py::def("collectElements", collectElements);
	py::def("collectCharges", collectCharges);
	py::def("collectMolarMasses", collectMolarMasses);

	export_std_vector<Species>("SpeciesVector");
}

} // namespace Reaktor
