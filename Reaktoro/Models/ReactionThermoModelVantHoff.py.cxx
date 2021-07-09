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
#include <Reaktoro/pybind11.hxx>

// Reaktoro includes
#include <Reaktoro/Models/ReactionThermoModelVantHoff.hpp>
using namespace Reaktoro;

void exportReactionThermoModelVantHoff(py::module& m)
{
    py::class_<ReactionThermoModelParamsVantHoff>(m, "ReactionThermoModelParamsVantHoff")
        .def_readwrite("lgKr", &ReactionThermoModelParamsVantHoff::lgKr)
        .def_readwrite("dHr", &ReactionThermoModelParamsVantHoff::dHr)
        .def_readwrite("Tr", &ReactionThermoModelParamsVantHoff::Tr)
        .def_readwrite("Pr", &ReactionThermoModelParamsVantHoff::Pr)
        ;

    m.def("ReactionThermoModelVantHoff", ReactionThermoModelVantHoff);
}
