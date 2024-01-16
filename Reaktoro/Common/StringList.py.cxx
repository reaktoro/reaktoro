// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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
#include <Reaktoro/Common/StringList.hpp>
using namespace Reaktoro;

void exportStringList(py::module& m)
{
    const auto __getitem__ = [](const StringList& self, std::size_t i) { return self[i]; };
    const auto __setitem__ = [](StringList& self, std::size_t i, std::string str) { self[i] = str; };

    py::class_<StringList>(m, "StringList")
        .def(py::init<>())
        .def(py::init<std::initializer_list<std::string>>())
        .def(py::init<std::vector<std::string>>())
        .def(py::init<const char*>())
        .def(py::init<const char*, char>())
        .def(py::init<std::string>())
        .def(py::init<std::string, char>())
        .def("empty", &StringList::empty)
        .def("size", &StringList::size)
        .def("data", &StringList::data, return_internal_ref)
        .def("__getitem__", __getitem__, return_internal_ref)
        .def("__setitem__", __setitem__)
        .def("__iter__", [](StringList const& self) { return py::make_iterator(self.begin(), self.end()); }, py::keep_alive<0, 1>()) // keep object alive while iterator exists;
        ;

    py::implicitly_convertible<std::string, StringList>();
    py::implicitly_convertible<std::vector<std::string>, StringList>();
}
