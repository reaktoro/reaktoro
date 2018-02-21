// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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

// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Core/Element.hpp>
#include <Reaktoro/Core/Species.hpp>

namespace Reaktoro {

void exportSpecies(py::module& m)
{
    py::class_<Species>(m, "Species")
        .def(py::init<>())
        .def("setName", &Species::setName)
        .def("setFormula", &Species::setFormula)
        .def("setElements", &Species::setElements)
        .def("numElements", &Species::numElements)
        .def("name", &Species::name)
        .def("formula", &Species::formula)
        .def("elements", &Species::elements, py::return_value_policy::reference_internal)
        .def("molarMass", &Species::molarMass)
        .def("charge", &Species::charge)
        .def("elementCoefficient", &Species::elementCoefficient)
        ;
}

} // namespace Reaktoro
