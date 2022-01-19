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
#include <Reaktoro/Thermodynamics/Aqueous/ActivityModelDrummond.hpp>
using namespace Reaktoro;

void exportActivityModelDrummond(py::module& m)
{
    py::class_<ActivityModelDrummondParams>(m, "ActivityModelDrummondParams")
        .def(py::init<>())
        .def_readwrite("a1", &ActivityModelDrummondParams::a1)
        .def_readwrite("a2", &ActivityModelDrummondParams::a2)
        .def_readwrite("a3", &ActivityModelDrummondParams::a3)
        .def_readwrite("a4", &ActivityModelDrummondParams::a4)
        .def_readwrite("a5", &ActivityModelDrummondParams::a5)
        ;

    m.def("ActivityModelDrummond", py::overload_cast<String>(ActivityModelDrummond));
    m.def("ActivityModelDrummond", py::overload_cast<String, ActivityModelDrummondParams>(ActivityModelDrummond));
}
