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
#include <Reaktoro/Thermodynamics/Surface/ActivityModelSurfaceComplexation.hpp>

using namespace Reaktoro;

void exportActivityModelSurfaceComplexation(py::module& m)
{
    m.def("ActivityModelSurfaceComplexationNoDDL", ActivityModelSurfaceComplexationNoDDL);
    m.def("ActivityModelSurfaceComplexationWithDDL", ActivityModelSurfaceComplexationWithDDL);
    m.def("ActivityModelDDL", ActivityModelDDL());
    m.def("ActivityModelDDL", py::overload_cast<ActivityModelDDLParams>(ActivityModelDDL));

    py::class_<ActivityModelSurfaceComplexationParams>(m, "ActivityModelSurfaceComplexationParams")
        .def(py::init<>())
        .def_readwrite("surface", &ActivityModelSurfaceComplexationParams::surface)
        ;

    py::class_<ActivityModelDDLParams>(m, "ActivityModelDDLParams")
        .def(py::init<>())
        .def_readwrite("enr", &ActivityModelDDLParams::enr)
        .def_readwrite("thickness", &ActivityModelDDLParams::thickness)
        ;

    py::class_<ActivityModelDDLDonnanParams>(m, "ActivityModelDDLDonnanParams")
        .def(py::init<>())
        .def_readwrite("ddl", &ActivityModelDDLDonnanParams::ddl)
        .def_readwrite("debye_lengths", &ActivityModelDDLDonnanParams::debye_lengths)
        .def_readwrite("limit", &ActivityModelDDLDonnanParams::limit)
        .def_readwrite("viscosity", &ActivityModelDDLDonnanParams::viscosity)
        ;
}