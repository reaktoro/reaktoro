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
#include <pybind11/eigen.h>
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Core/ActivityProps.hpp>
using namespace Reaktoro;

void exportActivityProps(py::module& m)
{
    py::class_<ActivityProps>(m, "ActivityProps")
        .def(py::init<>())
        .def_readwrite("Vex", &ActivityProps::Vex)
        .def_readwrite("VexT", &ActivityProps::VexT)
        .def_readwrite("VexP", &ActivityProps::VexP)
        .def_readwrite("Gex", &ActivityProps::Gex)
        .def_readwrite("Hex", &ActivityProps::Hex)
        .def_readwrite("Cpex", &ActivityProps::Cpex)
        .def_readwrite("Cvex", &ActivityProps::Cvex)
        .def_readwrite("ln_g", &ActivityProps::ln_g)
        .def_readwrite("ln_a", &ActivityProps::ln_a)
        ;

    py::class_<ActivityPropsConstRef>(m, "ActivityPropsConstRef")
        .def_property_readonly("Vex", [](const ActivityPropsConstRef& self) { return self.Vex; })
        .def_property_readonly("VexT", [](const ActivityPropsConstRef& self) { return self.VexT; })
        .def_property_readonly("VexP", [](const ActivityPropsConstRef& self) { return self.VexP; })
        .def_property_readonly("Gex", [](const ActivityPropsConstRef& self) { return self.Gex; })
        .def_property_readonly("Hex", [](const ActivityPropsConstRef& self) { return self.Hex; })
        .def_property_readonly("Cpex", [](const ActivityPropsConstRef& self) { return self.Cpex; })
        .def_property_readonly("Cvex", [](const ActivityPropsConstRef& self) { return self.Cvex; })
        .def_readonly("ln_g", &ActivityPropsConstRef::ln_g)
        .def_readonly("ln_a", &ActivityPropsConstRef::ln_a)
        ;

    py::class_<ActivityArgs>(m, "ActivityArgs")
        .def_property_readonly("T", [](const ActivityArgs& self) { return self.T; })
        .def_property_readonly("P", [](const ActivityArgs& self) { return self.P; })
        .def_property_readonly("x", [](const ActivityArgs& self) { return self.x; })
        .def_property_readonly("extra", [](const ActivityArgs& self) { return self.extra; })
        ;
}
