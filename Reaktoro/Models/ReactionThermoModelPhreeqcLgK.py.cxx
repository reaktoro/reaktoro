// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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
#include <pybind11/functional.h>
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Models/ReactionThermoModelPhreeqcLgK.hpp>
using namespace Reaktoro;

void exportReactionThermoModelPhreeqcLgK(py::module& m)
{
    py::class_<ReactionThermoModelParamsPhreeqcLgK>(m, "ReactionThermoModelParamsPhreeqcLgK")
        .def_readwrite("A1", &ReactionThermoModelParamsPhreeqcLgK::A1)
        .def_readwrite("A2", &ReactionThermoModelParamsPhreeqcLgK::A2)
        .def_readwrite("A3", &ReactionThermoModelParamsPhreeqcLgK::A3)
        .def_readwrite("A4", &ReactionThermoModelParamsPhreeqcLgK::A4)
        .def_readwrite("A5", &ReactionThermoModelParamsPhreeqcLgK::A5)
        .def_readwrite("A6", &ReactionThermoModelParamsPhreeqcLgK::A6)
        .def_readwrite("Pr", &ReactionThermoModelParamsPhreeqcLgK::Pr)
        ;

    m.def("ReactionThermoModelPhreeqcLgK", ReactionThermoModelPhreeqcLgK);
}
