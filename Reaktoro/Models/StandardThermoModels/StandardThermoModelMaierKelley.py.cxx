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
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelMaierKelley.hpp>
using namespace Reaktoro;

void exportStandardThermoModelMaierKelley(py::module& m)
{
    py::class_<StandardThermoModelParamsMaierKelley>(m, "StandardThermoModelParamsMaierKelley")
        .def(py::init<>())
        .def_readwrite("Gf",       &StandardThermoModelParamsMaierKelley::Gf)
        .def_readwrite("Hf",       &StandardThermoModelParamsMaierKelley::Hf)
        .def_readwrite("Sr",       &StandardThermoModelParamsMaierKelley::Sr)
        .def_readwrite("Vr",       &StandardThermoModelParamsMaierKelley::Vr)
        .def_readwrite("a",        &StandardThermoModelParamsMaierKelley::a)
        .def_readwrite("b",        &StandardThermoModelParamsMaierKelley::b)
        .def_readwrite("c",        &StandardThermoModelParamsMaierKelley::c)
        .def_readwrite("Tmax",     &StandardThermoModelParamsMaierKelley::Tmax)
        ;

    m.def("StandardThermoModelMaierKelley", StandardThermoModelMaierKelley);
}
