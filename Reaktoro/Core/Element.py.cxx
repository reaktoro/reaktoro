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
    auto createElement = [](
        String symbol,
        String name,
        Index atomic_number,
        double atomic_weight,
        double electronegativity,
        Strings tags)
    {
        return Element({symbol, name, atomic_number, atomic_weight, electronegativity, tags});
    };

    py::class_<Element>(m, "Element")
        .def(py::init<>())
        .def(py::init<String>())
        .def(py::init(createElement),
            py::arg("symbol"),
            py::arg("name"),
            py::arg("atomic_number"),
            py::arg("atomic_weight"),
            py::arg("electronegativity"),
            py::arg("tags"))
        .def("clone", &Element::clone)
        .def("withSymbol", &Element::withSymbol)
        .def("withName", &Element::withName)
        .def("withAtomicNumber", &Element::withAtomicNumber)
        .def("withAtomicWeight", &Element::withAtomicWeight)
        .def("withMolarMass", &Element::withMolarMass)
        .def("withElectronegativity", &Element::withElectronegativity)
        .def("withTags", &Element::withTags)
        .def("symbol", &Element::symbol)
        .def("name", &Element::name)
        .def("atomicNumber", &Element::atomicNumber)
        .def("atomicWeight", &Element::atomicWeight)
        .def("molarMass", &Element::molarMass)
        .def("electronegativity", &Element::electronegativity)
        .def("tags", &Element::tags, py::return_value_policy::reference_internal)
        ;
}
