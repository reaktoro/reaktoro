// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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
#include <Reaktoro/Core/ElementList.hpp>
using namespace Reaktoro;

void exportElementList(py::module& m)
{
    auto __getitem__ = [](const ElementList& self, Index i)
    {
        return self[i];
    };

    py::class_<ElementList>(m, "ElementList")
        .def(py::init<>())
        .def(py::init<const Vec<Element>&>())
        .def("append", &ElementList::append)
        .def("data", &ElementList::data)
        .def("size", &ElementList::size)
        .def("find", &ElementList::find)
        .def("findWithSymbol", &ElementList::findWithSymbol)
        .def("findWithName", &ElementList::findWithName)
        .def("index", &ElementList::index)
        .def("indexWithSymbol", &ElementList::indexWithSymbol)
        .def("indexWithName", &ElementList::indexWithName)
        .def("get", &ElementList::get, py::return_value_policy::reference_internal)
        .def("getWithSymbol", &ElementList::getWithSymbol, py::return_value_policy::reference_internal)
        .def("getWithName", &ElementList::getWithName, py::return_value_policy::reference_internal)
        .def("withSymbols", &ElementList::withSymbols)
        .def("withNames", &ElementList::withNames)
        .def("withTag", &ElementList::withTag)
        .def("withoutTag", &ElementList::withoutTag)
        .def("withTags", &ElementList::withTags)
        .def("withoutTags", &ElementList::withoutTags)
        .def("__getitem__", __getitem__, py::return_value_policy::reference_internal)
        .def("__iter__", [](const ElementList& self) { return py::make_iterator(self.begin(), self.end()); }, py::keep_alive<0, 1>()); // keep object alive while iterator exists;
        ;
}
