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

        .def("value", [](Param& self, const real& val) { return self.value(val); }, return_internal_ref)
        .def("value", [](const Param& self) { return self.value(); }, return_internal_ref)

        .def("lowerbound", [](Param& self, const real& val) { return self.lowerbound(val); }, return_internal_ref)
        .def("lowerbound", [](const Param& self) { return self.lowerbound(); })

        .def("upperbound", [](Param& self, const real& val) { return self.upperbound(val); }, return_internal_ref)
        .def("upperbound", [](const Param& self) { return self.upperbound(); })

        .def("isconst", [](Param& self, bool val) { return self.isconst(val); }, return_internal_ref)
        .def("isconst", [](const Param& self) { return self.isconst(); })

        .def_static("Constant", &Param::Constant)
        ;

    py::implicitly_convertible<real, Param>();
    py::implicitly_convertible<double, Param>();
}
