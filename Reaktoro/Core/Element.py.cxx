// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Core/Element.hpp>
using namespace Reaktoro;

void exportElement(py::module& m)
{
    auto createElement = [](String symbol, double molar_mass, String name, Strings tags)
    {
        return Element({symbol, molar_mass, name, tags});
    };

    py::class_<Element>(m, "Element")
        .def(py::init<>())
        .def(py::init<String>())
        .def(py::init(createElement),
            py::arg("symbol"),
            py::arg("molar_mass"),
            py::arg("name"),
            py::arg("tags"))
        .def("clone", &Element::clone)
        .def("withSymbol", &Element::withSymbol)
        .def("withMolarMass", &Element::withMolarMass)
        .def("withName", &Element::withName)
        .def("withTags", &Element::withTags)
        .def("symbol", &Element::symbol)
        .def("molarMass", &Element::molarMass)
        .def("name", &Element::name)
        .def("tags", &Element::tags, py::return_value_policy::reference_internal)
        ;
}
