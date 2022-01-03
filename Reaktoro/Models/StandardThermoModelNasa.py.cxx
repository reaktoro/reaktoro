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
#include <Reaktoro/Models/StandardThermoModelNasa.hpp>
using namespace Reaktoro;

void exportStandardThermoModelNasa(py::module& m)
{
    py::class_<StandardThermoModelParamsNasa::TemperatureInterval>(m, "StandardThermoModelParamsNasaTemperatureInterval")
        .def_readwrite("Tmin",           &StandardThermoModelParamsNasa::TemperatureInterval::Tmin)
        .def_readwrite("Tmax",           &StandardThermoModelParamsNasa::TemperatureInterval::Tmax)
        .def_readwrite("label",          &StandardThermoModelParamsNasa::TemperatureInterval::label)
        .def_readwrite("aggregatestate", &StandardThermoModelParamsNasa::TemperatureInterval::aggregatestate)
        .def_readwrite("a1",             &StandardThermoModelParamsNasa::TemperatureInterval::a1)
        .def_readwrite("a2",             &StandardThermoModelParamsNasa::TemperatureInterval::a2)
        .def_readwrite("a3",             &StandardThermoModelParamsNasa::TemperatureInterval::a3)
        .def_readwrite("a4",             &StandardThermoModelParamsNasa::TemperatureInterval::a4)
        .def_readwrite("a5",             &StandardThermoModelParamsNasa::TemperatureInterval::a5)
        .def_readwrite("a6",             &StandardThermoModelParamsNasa::TemperatureInterval::a6)
        .def_readwrite("a7",             &StandardThermoModelParamsNasa::TemperatureInterval::a7)
        .def_readwrite("b1",             &StandardThermoModelParamsNasa::TemperatureInterval::b1)
        .def_readwrite("b2",             &StandardThermoModelParamsNasa::TemperatureInterval::b2)
        ;

    py::class_<StandardThermoModelParamsNasa>(m, "StandardThermoModelParamsNasa")
        .def_readwrite("dHf",       &StandardThermoModelParamsNasa::dHf)
        .def_readwrite("dH0",       &StandardThermoModelParamsNasa::dH0)
        .def_readwrite("H0",        &StandardThermoModelParamsNasa::H0)
        .def_readwrite("T0",        &StandardThermoModelParamsNasa::T0)
        .def_readwrite("intervals", &StandardThermoModelParamsNasa::intervals)
        ;

    m.def("StandardThermoModelNasa", StandardThermoModelNasa);
}
