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
#include <Reaktoro/pybind11.hxx>

// Reaktoro includes
#include <Reaktoro/Core/ElementList.hpp>
using namespace Reaktoro;

template<typename ElementListType>
auto createCommonElementList(py::module& m, const char* name)
{
    return
    py::class_<ElementListType>(m, name)
        .def("data", &ElementListType::data)
        .def("size", &ElementListType::size)
        .def("find", &ElementListType::find)
        .def("findWithSymbol", &ElementListType::findWithSymbol)
        .def("findWithName", &ElementListType::findWithName)
        .def("index", &ElementListType::index)
        .def("indexWithSymbol", &ElementListType::indexWithSymbol)
        .def("indexWithName", &ElementListType::indexWithName)
        .def("get", &ElementListType::get, py::return_value_policy::reference_internal)
        .def("getWithSymbol", &ElementListType::getWithSymbol, py::return_value_policy::reference_internal)
        .def("getWithName", &ElementListType::getWithName, py::return_value_policy::reference_internal)
        .def("withSymbols", &ElementListType::withSymbols)
        .def("withNames", &ElementListType::withNames)
        .def("withTag", &ElementListType::withTag)
        .def("withoutTag", &ElementListType::withoutTag)
        .def("withTags", &ElementListType::withTags)
        .def("withoutTags", &ElementListType::withoutTags)
        .def("__getitem__", [](const ElementListType& self, Index i) { return self[i]; }, py::return_value_policy::reference_internal)
        .def("__iter__", [](const ElementListType& self) { return py::make_iterator(self.begin(), self.end()); }, py::keep_alive<0, 1>()); // keep object alive while iterator exists;
        ;
}

void exportElementList(py::module& m)
{
    createCommonElementList<ElementList>(m, "ElementList")
        .def(py::init<>())
        .def(py::init<const Vec<Element>&>())
        .def("append", &ElementList::append)
        ;

    createCommonElementList<ElementListConstRef>(m, "ElementListConstRef");
}
