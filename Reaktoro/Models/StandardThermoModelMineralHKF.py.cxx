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
#include <Reaktoro/Models/StandardThermoModelMineralHKF.hpp>
using namespace Reaktoro;

void exportStandardThermoModelMineralHKF(py::module& m)
{
    py::class_<StandardThermoModelParamsMineralHKF>(m, "StandardThermoModelParamsMineralHKF")
        .def_readwrite("Gf",     &StandardThermoModelParamsMineralHKF::Gf)
        .def_readwrite("Hf",     &StandardThermoModelParamsMineralHKF::Hf)
        .def_readwrite("Sr",     &StandardThermoModelParamsMineralHKF::Sr)
        .def_readwrite("Vr",     &StandardThermoModelParamsMineralHKF::Vr)
        .def_readwrite("ntr",    &StandardThermoModelParamsMineralHKF::ntr)
        .def_readwrite("a",      &StandardThermoModelParamsMineralHKF::a)
        .def_readwrite("b",      &StandardThermoModelParamsMineralHKF::b)
        .def_readwrite("c",      &StandardThermoModelParamsMineralHKF::c)
        .def_readwrite("Ttr",    &StandardThermoModelParamsMineralHKF::Ttr)
        .def_readwrite("Htr",    &StandardThermoModelParamsMineralHKF::Htr)
        .def_readwrite("Vtr",    &StandardThermoModelParamsMineralHKF::Vtr)
        .def_readwrite("dPdTtr", &StandardThermoModelParamsMineralHKF::dPdTtr)
        .def_readwrite("Tmax",   &StandardThermoModelParamsMineralHKF::Tmax)
        ;

    m.def("StandardThermoModelMineralHKF", StandardThermoModelMineralHKF);
}
