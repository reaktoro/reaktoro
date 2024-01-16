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
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelInterpolation.hpp>
using namespace Reaktoro;

void exportStandardThermoModelInterpolation(py::module& m)
{
    py::class_<StandardThermoModelParamsInterpolation>(m, "StandardThermoModelParamsInterpolation")
        .def(py::init<>())
        .def_readwrite("temperatures", &StandardThermoModelParamsInterpolation::temperatures)
        .def_readwrite("pressures",    &StandardThermoModelParamsInterpolation::pressures)
        .def_readwrite("G0",           &StandardThermoModelParamsInterpolation::G0)
        .def_readwrite("H0",           &StandardThermoModelParamsInterpolation::H0)
        .def_readwrite("V0",           &StandardThermoModelParamsInterpolation::V0)
        .def_readwrite("VT0",          &StandardThermoModelParamsInterpolation::VT0)
        .def_readwrite("VP0",          &StandardThermoModelParamsInterpolation::VP0)
        .def_readwrite("Cp0",          &StandardThermoModelParamsInterpolation::Cp0)
        ;

    m.def("StandardThermoModelInterpolation", StandardThermoModelInterpolation);
}
