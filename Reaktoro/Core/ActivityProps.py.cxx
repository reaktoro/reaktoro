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
#include <Reaktoro/Core/ActivityProps.hpp>
using namespace Reaktoro;

void exportActivityProps(py::module& m)
{
    py::class_<ActivityProps>(m, "ActivityProps")
        .def(py::init<>())
        .def_readwrite("Vx", &ActivityProps::Vx)
        .def_readwrite("VxT", &ActivityProps::VxT)
        .def_readwrite("VxP", &ActivityProps::VxP)
        .def_readwrite("Gx", &ActivityProps::Gx)
        .def_readwrite("Hx", &ActivityProps::Hx)
        .def_readwrite("Cpx", &ActivityProps::Cpx)
        .def_readwrite("ln_g", &ActivityProps::ln_g)
        .def_readwrite("ln_a", &ActivityProps::ln_a)
        .def_readwrite("extra", &ActivityProps::extra)
        .def("create", &ActivityProps::create)
        ;

    #define get(field) [](const ActivityPropsRef& self) { return self.field; }
    #define set(field) [](ActivityPropsRef& self, decltype(ActivityPropsRef::field)& val) { self.field = val; }

    py::class_<ActivityPropsRef>(m, "ActivityPropsRef")
        .def_property("Vx", get(Vx), set(Vx))
        .def_property("VxT", get(VxT), set(VxT))
        .def_property("VxP", get(VxP), set(VxP))
        .def_property("Gx", get(Gx), set(Gx))
        .def_property("Hx", get(Hx), set(Hx))
        .def_property("Cpx", get(Cpx), set(Cpx))
        .def_readwrite("ln_g", &ActivityPropsRef::ln_g)
        .def_readwrite("ln_a", &ActivityPropsRef::ln_a)
        .def_property("extra", get(extra), set(extra))
        ;

    #undef get
    #undef set

    #define get(field) [](const ActivityPropsConstRef& self) { return self.field; }

    py::class_<ActivityPropsConstRef>(m, "ActivityPropsConstRef")
        .def_property_readonly("Vx", get(Vx))
        .def_property_readonly("VxT", get(VxT))
        .def_property_readonly("VxP", get(VxP))
        .def_property_readonly("Gx", get(Gx))
        .def_property_readonly("Hx", get(Hx))
        .def_property_readonly("Cpx", get(Cpx))
        .def_readonly("ln_g", &ActivityPropsConstRef::ln_g)
        .def_readonly("ln_a", &ActivityPropsConstRef::ln_a)
        .def_property_readonly("extra", get(extra))
        ;
}
