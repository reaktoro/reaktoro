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
#include <Reaktoro/Models/StandardThermoModels/ReactionStandardThermoModelPhreeqcLgK.hpp>
using namespace Reaktoro;

void exportReactionStandardThermoModelPhreeqcLgK(py::module& m)
{
    py::class_<ReactionStandardThermoModelParamsPhreeqcLgK>(m, "ReactionStandardThermoModelParamsPhreeqcLgK")
        .def(py::init<>())
        .def_readwrite("A1", &ReactionStandardThermoModelParamsPhreeqcLgK::A1)
        .def_readwrite("A2", &ReactionStandardThermoModelParamsPhreeqcLgK::A2)
        .def_readwrite("A3", &ReactionStandardThermoModelParamsPhreeqcLgK::A3)
        .def_readwrite("A4", &ReactionStandardThermoModelParamsPhreeqcLgK::A4)
        .def_readwrite("A5", &ReactionStandardThermoModelParamsPhreeqcLgK::A5)
        .def_readwrite("A6", &ReactionStandardThermoModelParamsPhreeqcLgK::A6)
        .def_readwrite("Pr", &ReactionStandardThermoModelParamsPhreeqcLgK::Pr)
        ;

    m.def("ReactionStandardThermoModelPhreeqcLgK", ReactionStandardThermoModelPhreeqcLgK);
}
