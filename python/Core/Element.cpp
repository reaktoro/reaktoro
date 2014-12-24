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

#include "Element.hpp"

// Boost includes
#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

using namespace boost::python;

// Reaktor includes
#include <Reaktor/Core/Element.hpp>

// PyReator includes
#include <python/Utils/Converters.hpp>

namespace Reaktor {

auto createElement(std::string name, double molar_mass) -> boost::shared_ptr<Element>
{
	ElementData data;
	data.name = name;
	data.molar_mass = molar_mass;
	return boost::make_shared<Element>(data);
}

auto export_Element() -> void
{
	class_<ElementData>("ElementData")
		.def_readwrite("name", &ElementData::name)
		.def_readwrite("molar_mass", &ElementData::molar_mass)
		;

	class_<Element>("Element")
		.def(init<>())
		.def(init<const ElementData&>())
		.def("__init__", make_constructor(createElement, default_call_policies(), (arg("name"), arg("molar_mass"))))
		.def("name", &Element::name)
		.def("molarMass", &Element::molarMass)
		;

	export_std_vector<Element>("ElementVector");
}

} // namespace Reaktor
