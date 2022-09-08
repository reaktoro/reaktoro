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
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelHKF.hpp>
using namespace Reaktoro;

void exportStandardThermoModelHKF(py::module& m)
{
    py::class_<StandardThermoModelParamsHKF>(m, "StandardThermoModelParamsHKF")
        .def_readwrite("Gf",     &StandardThermoModelParamsHKF::Gf)
        .def_readwrite("Hf",     &StandardThermoModelParamsHKF::Hf)
        .def_readwrite("Sr",     &StandardThermoModelParamsHKF::Sr)
        .def_readwrite("a1",     &StandardThermoModelParamsHKF::a1)
        .def_readwrite("a2",     &StandardThermoModelParamsHKF::a2)
        .def_readwrite("a3",     &StandardThermoModelParamsHKF::a3)
        .def_readwrite("a4",     &StandardThermoModelParamsHKF::a4)
        .def_readwrite("c1",     &StandardThermoModelParamsHKF::c1)
        .def_readwrite("c2",     &StandardThermoModelParamsHKF::c2)
        .def_readwrite("wref",   &StandardThermoModelParamsHKF::wref)
        .def_readwrite("charge", &StandardThermoModelParamsHKF::charge)
        .def_readwrite("Tmax",   &StandardThermoModelParamsHKF::Tmax)
        ;

    m.def("StandardThermoModelHKF", StandardThermoModelHKF);
}
