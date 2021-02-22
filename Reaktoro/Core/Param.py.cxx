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
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Core/Param.hpp>
using namespace Reaktoro;

void exportParam(py::module& m)
{
    const auto return_internal_ref = py::return_value_policy::reference_internal;

    py::class_<Param>(m, "Param")
        .def(py::init<>())
        .def(py::init<const real&>())
        .def(py::init<const double&>())
        .def(py::init<const long&>())
        .def(py::init<const Param&>())
        .def(py::init<const String&, const real&>())
        .def(py::init<const String&, double>())
        .def(py::init<const String&, long>())

        .def("assign", &Param::assign, return_internal_ref);

        .def("value", py::overload_cast<const real&>(&Param::value), return_internal_ref)
        .def("value", py::overload_cast<>(&Param::value, py::const_), return_internal_ref)

        .def("id", py::overload_cast<String>(&Param::id), return_internal_ref)
        .def("id", py::overload_cast<>(&Param::id, py::const_), return_internal_ref)

        .def("lowerbound", py::overload_cast<double>(&Param::lowerbound), return_internal_ref)
        .def("lowerbound", py::overload_cast<>(&Param::lowerbound, py::const_))

        .def("upperbound", py::overload_cast<double>(&Param::upperbound), return_internal_ref)
        .def("upperbound", py::overload_cast<>(&Param::upperbound, py::const_))

        .def("isconst", py::overload_cast<bool>(&Param::isconst), return_internal_ref)
        .def("isconst", py::overload_cast<>(&Param::isconst, py::const_))

        .def_static("Constant", &Param::Constant)
        ;

    py::implicitly_convertible<real, Param>();
    py::implicitly_convertible<double, Param>();
}
