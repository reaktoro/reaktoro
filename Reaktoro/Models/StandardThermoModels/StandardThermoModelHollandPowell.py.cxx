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
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelHollandPowell.hpp>
using namespace Reaktoro;

void exportStandardThermoModelHollandPowell(py::module& m)
{
    py::class_<StandardThermoModelParamsHollandPowell>(m, "StandardThermoModelParamsHollandPowell")
        .def_readwrite("Gf",       &StandardThermoModelParamsHollandPowell::Gf)
        .def_readwrite("Hf",       &StandardThermoModelParamsHollandPowell::Hf)
        .def_readwrite("Sr",       &StandardThermoModelParamsHollandPowell::Sr)
        .def_readwrite("Vr",       &StandardThermoModelParamsHollandPowell::Vr)
        .def_readwrite("a",        &StandardThermoModelParamsHollandPowell::a)
        .def_readwrite("b",        &StandardThermoModelParamsHollandPowell::b)
        .def_readwrite("c",        &StandardThermoModelParamsHollandPowell::c)
        .def_readwrite("d",        &StandardThermoModelParamsHollandPowell::d)
        .def_readwrite("alpha0",   &StandardThermoModelParamsHollandPowell::alpha0)
        .def_readwrite("kappa0",   &StandardThermoModelParamsHollandPowell::kappa0)
        .def_readwrite("kappa0p",  &StandardThermoModelParamsHollandPowell::kappa0p)
        .def_readwrite("kappa0pp", &StandardThermoModelParamsHollandPowell::kappa0pp)
        .def_readwrite("numatoms", &StandardThermoModelParamsHollandPowell::numatoms)
        .def_readwrite("Tmax",     &StandardThermoModelParamsHollandPowell::Tmax)
        ;

    m.def("StandardThermoModelHollandPowell", StandardThermoModelHollandPowell);
}
