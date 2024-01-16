// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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
#include <Reaktoro/Models/ActivityModels/ActivityModelDebyeHuckel.hpp>
using namespace Reaktoro;

void exportActivityModelDebyeHuckel(py::module& m)
{
    py::class_<ActivityModelDebyeHuckelParams>(m, "ActivityModelDebyeHuckelParams")
        .def(py::init<>())
        .def_readwrite("aiondefault", &ActivityModelDebyeHuckelParams::aiondefault)
        .def_readwrite("biondefault", &ActivityModelDebyeHuckelParams::biondefault)
        .def_readwrite("bneutraldefault", &ActivityModelDebyeHuckelParams::bneutraldefault)
        .def_readwrite("aions", &ActivityModelDebyeHuckelParams::aions)
        .def_readwrite("bions", &ActivityModelDebyeHuckelParams::bions)
        .def_readwrite("bneutrals", &ActivityModelDebyeHuckelParams::bneutrals)
        .def("aion", &ActivityModelDebyeHuckelParams::aion)
        .def("bion", &ActivityModelDebyeHuckelParams::bion)
        .def("bneutral", &ActivityModelDebyeHuckelParams::bneutral)
        .def("setLimitingLaw", &ActivityModelDebyeHuckelParams::setLimitingLaw)
        .def("setKielland", &ActivityModelDebyeHuckelParams::setKielland)
        .def("setWATEQ4F", &ActivityModelDebyeHuckelParams::setWATEQ4F)
        .def("setPHREEQC", &ActivityModelDebyeHuckelParams::setPHREEQC)
        ;

    m.def("ActivityModelDebyeHuckel", py::overload_cast<>(ActivityModelDebyeHuckel));
    m.def("ActivityModelDebyeHuckel", py::overload_cast<ActivityModelDebyeHuckelParams>(ActivityModelDebyeHuckel));
    m.def("ActivityModelDebyeHuckelLimitingLaw", ActivityModelDebyeHuckelLimitingLaw);
    m.def("ActivityModelDebyeHuckelKielland", ActivityModelDebyeHuckelKielland);
    m.def("ActivityModelDebyeHuckelPHREEQC", ActivityModelDebyeHuckelPHREEQC);
    m.def("ActivityModelDebyeHuckelWATEQ4F", ActivityModelDebyeHuckelWATEQ4F);
}
