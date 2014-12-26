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

#include "PyElement.hpp"

// Boost includes
#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
namespace py = boost::python;

// Reaktor includes
#include <Reaktor/Reaktor.hpp>

// PyReator includes
#include <PyReaktor/Utils/PyConverters.hpp>

namespace Reaktor {

auto createElement(std::string name, double molar_mass) -> boost::shared_ptr<Element>
{
	ElementData data;
	data.name = std::move(name);
	data.molar_mass = std::move(molar_mass);
	return boost::make_shared<Element>(data);
}

auto export_Element() -> void
{
	py::class_<ElementData>("ElementData")
		.def_readwrite("name", &ElementData::name)
		.def_readwrite("molar_mass", &ElementData::molar_mass)
		;

	py::class_<Element>("Element")
		.def(py::init<>())
		.def(py::init<const ElementData&>())
		.def("__init__", py::make_constructor(createElement, py::default_call_policies(),
			(py::arg("name"), py::arg("molar_mass"))))
		.def("name", &Element::name)
		.def("molarMass", &Element::molarMass)
		;

	export_std_vector<Element>("ElementVector");
}

} // namespace Reaktor
