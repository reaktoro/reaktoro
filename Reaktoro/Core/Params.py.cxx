// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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
    const auto return_internal_ref = py::return_value_policy::reference_internal;

    py::class_<Params>(m, "Params")
        .def(py::init<>())
        .def("at", &Params::at, return_internal_ref)
        .def("get", &Params::get, return_internal_ref)
        .def("set", py::overload_cast<const String&, const Params&>(&Params::set))
        .def("set", py::overload_cast<const String&, const Param&>(&Params::set))
        ;
}
