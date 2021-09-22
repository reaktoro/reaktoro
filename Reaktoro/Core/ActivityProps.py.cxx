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
        .def_readwrite("extra", &ActivityProps::extra)
        .def("create", &ActivityProps::create)
        ;

    #define get(field) [](const ActivityPropsRef& self) { return self.field; }
    #define set(field) [](ActivityPropsRef& self, decltype(ActivityPropsRef::field)& val) { self.field = val; }

    py::class_<ActivityPropsRef>(m, "ActivityPropsRef")
        .def_property("Vex", get(Vex), set(Vex))
        .def_property("VexT", get(VexT), set(VexT))
        .def_property("VexP", get(VexP), set(VexP))
        .def_property("Gex", get(Gex), set(Gex))
        .def_property("Hex", get(Hex), set(Hex))
        .def_property("Cpex", get(Cpex), set(Cpex))
        .def_property("Cvex", get(Cvex), set(Cvex))
        .def_readwrite("ln_g", &ActivityPropsRef::ln_g)
        .def_readwrite("ln_a", &ActivityPropsRef::ln_a)
        .def_property("extra", get(extra), set(extra))
        ;

    #undef get
    #undef set

    #define get(field) [](const ActivityPropsConstRef& self) { return self.field; }

    py::class_<ActivityPropsConstRef>(m, "ActivityPropsConstRef")
        .def_property_readonly("Vex", get(Vex))
        .def_property_readonly("VexT", get(VexT))
        .def_property_readonly("VexP", get(VexP))
        .def_property_readonly("Gex", get(Gex))
        .def_property_readonly("Hex", get(Hex))
        .def_property_readonly("Cpex", get(Cpex))
        .def_property_readonly("Cvex", get(Cvex))
        .def_readonly("ln_g", &ActivityPropsConstRef::ln_g)
        .def_readonly("ln_a", &ActivityPropsConstRef::ln_a)
        .def_property_readonly("extra", get(extra))
        ;
}
