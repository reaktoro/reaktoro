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
#include <Reaktoro/Models/ActivityModels/ActivityModelExtendedUNIQUAC.hpp>
using namespace Reaktoro;

void exportActivityModelExtendedUNIQUAC(py::module& m)
{
    py::class_<ActivityModelParamsExtendedUNIQUAC>(m, "ActivityModelParamsExtendedUNIQUAC")
        .def(py::init<>())
        .def_readwrite("r", &ActivityModelParamsExtendedUNIQUAC::r, "The volume parameters :math:`r_i` for the aqueous species.")
        .def_readwrite("q", &ActivityModelParamsExtendedUNIQUAC::q, "The surface parameters :math:`r_i` for the aqueous species.")
        .def_readwrite("u", &ActivityModelParamsExtendedUNIQUAC::u, "The binary interaction parameters :math:`u_{ij}` among the aqueous species.")
        ;

    m.def("ActivityModelExtendedUNIQUAC", py::overload_cast<>(ActivityModelExtendedUNIQUAC));
    m.def("ActivityModelExtendedUNIQUAC", py::overload_cast<ActivityModelParamsExtendedUNIQUAC const&>(ActivityModelExtendedUNIQUAC));
    m.def("ActivityModelExtendedUNIQUAC", py::overload_cast<Params const&>(ActivityModelExtendedUNIQUAC));
}
