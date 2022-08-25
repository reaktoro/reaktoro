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
#include <Reaktoro/Thermodynamics/Surface/ActivityModelSorption.hpp>

using namespace Reaktoro;

void exportActivityModelSorption(py::module& m)
{
    m.def("ActivityModelSorptionNoDDL", ActivityModelSorptionNoDDL);
//    m.def("ActivityModelSorptionWithDDL", ActivityModelSorptionWithDDL);
//    m.def("ActivityModelSorptionWithEDL", ActivityModelSorptionWithEDL);
//    m.def("ActivityModelEDL", py::overload_cast<ActivityModelDDLParams>(ActivityModelEDL));

    py::class_<ActivityModelSorptionParams>(m, "ActivityModelSorptionParams")
        .def(py::init<>())
        .def_readwrite("surface", &ActivityModelSorptionParams::surface)
        .def_readwrite("site_tag", &ActivityModelSorptionParams::site_tag)
        .def_readwrite("output", &ActivityModelSorptionParams::output)
        ;

    py::class_<ActivityModelDDLParams>(m, "ActivityModelDDLParams")
        .def(py::init<>())
        .def_readwrite("enr", &ActivityModelDDLParams::enr)
        .def_readwrite("thickness", &ActivityModelDDLParams::thickness)
        ;

//    py::class_<ActivityModelDDLDonnanParams>(m, "ActivityModelDDLDonnanParams")
//        .def(py::init<>())
//        .def_readwrite("ddl", &ActivityModelDDLDonnanParams::ddl)
//        .def_readwrite("debye_lengths", &ActivityModelDDLDonnanParams::debye_lengths)
//        .def_readwrite("limit", &ActivityModelDDLDonnanParams::limit)
//        .def_readwrite("viscosity", &ActivityModelDDLDonnanParams::viscosity)
//        ;
}