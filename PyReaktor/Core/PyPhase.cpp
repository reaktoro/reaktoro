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

#include "PyPhase.hpp"

// Boost includes
#include <boost/python.hpp>
#include <boost/smart_ptr.hpp>
namespace py = boost::python;

// Reaktor includes
#include <Reaktor/Core/Species.hpp>
#include <Reaktor/Core/Phase.hpp>

// PyReator includes
#include <PyReaktor/Utils/PyConverters.hpp>

namespace Reaktor {

auto createPhase(std::string name, std::vector<Species> species) -> boost::shared_ptr<Phase>
{
	PhaseData data;
	data.name = std::move(name);
	data.species = std::move(species);
	return boost::make_shared<Phase>(data);
}

auto export_Phase() -> void
{
	py::class_<PhaseData>("PhaseData")
		.def_readwrite("name", &PhaseData::name)
		.def_readwrite("species", &PhaseData::species)
		;

	py::class_<Phase>("Phase")
		.def(py::init<>())
		.def(py::init<const PhaseData&>())
		.def("__init__", py::make_constructor(createPhase, py::default_call_policies(),
			(py::arg("name"), py::arg("species"))))
		.def("name", &Phase::name, py::return_value_policy<py::copy_const_reference>())
		.def("species", &Phase::species, py::return_value_policy<py::copy_const_reference>())
		;

	export_std_vector<Phase>("PhaseVector");
}

} // namespace Reaktor
