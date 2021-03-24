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
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Core/Params.hpp>
using namespace Reaktoro;

void exportParams(py::module& m)
{
    const auto return_internal_ref = py::return_value_policy::reference_internal;

    py::class_<Params>(m, "Params")
        .def(py::init<>())
        .def(py::init<const std::initializer_list<Param>&>())
        .def("clone", &Params::clone)
        .def("append", py::overload_cast<const Param&>(&Params::append), return_internal_ref)
        .def("append", py::overload_cast<const String&, const real&>(&Params::append), return_internal_ref)
        .def("resize", &Params::resize)
        .def("size", &Params::size)
        .def("find", &Params::find)
        .def("index", &Params::index)
        .def("get", py::overload_cast<const String&>(&Params::get), return_internal_ref)
        .def("get", py::overload_cast<const String&>(&Params::get, py::const_), return_internal_ref)
        .def("exists", &Params::exists)
        .def("data", &Params::data, return_internal_ref)
        .def("__getitem__", [](const Params& self, Index i) { return self[i]; })
        .def("__setitem__", [](Params& self, Index i, const Param& val) { self[i] = val; })
        .def("__getitem__", [](const Params& self, const String& id) { return self[id]; })
        .def("__setitem__", [](Params& self, const String& id, const Param& val) { self[id] = val; })
        .def("__iter__", [](const Params& self) { return py::make_iterator(self.begin(), self.end()); }, py::keep_alive<0, 1>()); // keep object alive while iterator exists;
        ;
}
