// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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
#include <Reaktoro/Core/Params.hpp>
using namespace Reaktoro;

void exportParams(py::module& m)
{
    // py::class_<Vec<Param>>(m, "Vec<Param>")
    //     .def(py::init<>())
    //     .def(py::init<const std::initializer_list<Param>&>())
    //     .def("clone", &Vec<Param>::clone)
    //     .def("append", py::overload_cast<const Param&>(&Vec<Param>::append), return_internal_ref)
    //     .def("append", py::overload_cast<const String&, const real&>(&Vec<Param>::append), return_internal_ref)
    //     .def("resize", &Vec<Param>::resize)
    //     .def("size", &Vec<Param>::size)
    //     .def("find", &Vec<Param>::find)
    //     .def("index", &Vec<Param>::index)
    //     .def("get", py::overload_cast<const String&>(&Vec<Param>::get), return_internal_ref)
    //     .def("get", py::overload_cast<const String&>(&Vec<Param>::get, py::const_), return_internal_ref)
    //     .def("exists", &Vec<Param>::exists)
    //     .def("data", &Vec<Param>::data, return_internal_ref)
    //     .def("__getitem__", [](const Vec<Param>& self, Index i) { return self[i]; })
    //     .def("__setitem__", [](Vec<Param>& self, Index i, const Param& val) { self[i] = val; })
    //     .def("__getitem__", [](const Vec<Param>& self, const String& id) { return self[id]; })
    //     .def("__setitem__", [](Vec<Param>& self, const String& id, const Param& val) { self[id] = val; })
    //     .def("__iter__", [](const Vec<Param>& self) { return py::make_iterator(self.begin(), self.end()); }, py::keep_alive<0, 1>()); // keep object alive while iterator exists;
    //     ;
}
