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

#include "PyEquilibriumProblem.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktor includes
#include <Reaktor/Reaktor.hpp>

namespace Reaktor {

auto export_EquilibriumProblem() -> void
{
	py::class_<EquilibriumProblem>("EquilibriumProblem", py::no_init)
        .def(py::init<const ChemicalSystem&>())
        .def(py::init<const ChemicalSystem&, const Partition&>())
        .def(py::init<const EquilibriumProblem&>())
        .def("setTemperature", &EquilibriumProblem::setTemperature, py::return_internal_reference<>())
        .def("setPressure", &EquilibriumProblem::setPressure, py::return_internal_reference<>())
        .def("setCharge", &EquilibriumProblem::setCharge, py::return_internal_reference<>())
        .def("setElementAmounts", &EquilibriumProblem::setElementAmounts, py::return_internal_reference<>())
        .def("temperature", &EquilibriumProblem::temperature)
        .def("pressure", &EquilibriumProblem::pressure)
        .def("charge", &EquilibriumProblem::charge)
        .def("elementAmounts", &EquilibriumProblem::elementAmounts, py::return_value_policy<py::copy_const_reference>())
        .def("componentAmounts", &EquilibriumProblem::componentAmounts)
        .def("balanceMatrix", &EquilibriumProblem::balanceMatrix, py::return_value_policy<py::copy_const_reference>())
        .def("components", &EquilibriumProblem::components, py::return_value_policy<py::copy_const_reference>())
        .def("system", &EquilibriumProblem::system, py::return_value_policy<py::copy_const_reference>())
        .def("partition", &EquilibriumProblem::partition, py::return_value_policy<py::copy_const_reference>())
		;

	py::implicitly_convertible<EquilibriumProblem, OptimumProblem>();
}

} // namespace Reaktor
