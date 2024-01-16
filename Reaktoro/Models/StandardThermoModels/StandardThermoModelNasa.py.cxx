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
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelNasa.hpp>
using namespace Reaktoro;

void exportStandardThermoModelNasa(py::module& m)
{
    py::class_<StandardThermoModelParamsNasa::Polynomial>(m, "StandardThermoModelParamsNasaPolynomial")
        .def(py::init<>())
        .def_readwrite("Tmin",  &StandardThermoModelParamsNasa::Polynomial::Tmin)
        .def_readwrite("Tmax",  &StandardThermoModelParamsNasa::Polynomial::Tmax)
        .def_readwrite("label", &StandardThermoModelParamsNasa::Polynomial::label)
        .def_readwrite("state", &StandardThermoModelParamsNasa::Polynomial::state)
        .def_readwrite("a1",    &StandardThermoModelParamsNasa::Polynomial::a1)
        .def_readwrite("a2",    &StandardThermoModelParamsNasa::Polynomial::a2)
        .def_readwrite("a3",    &StandardThermoModelParamsNasa::Polynomial::a3)
        .def_readwrite("a4",    &StandardThermoModelParamsNasa::Polynomial::a4)
        .def_readwrite("a5",    &StandardThermoModelParamsNasa::Polynomial::a5)
        .def_readwrite("a6",    &StandardThermoModelParamsNasa::Polynomial::a6)
        .def_readwrite("a7",    &StandardThermoModelParamsNasa::Polynomial::a7)
        .def_readwrite("b1",    &StandardThermoModelParamsNasa::Polynomial::b1)
        .def_readwrite("b2",    &StandardThermoModelParamsNasa::Polynomial::b2)
        ;

    py::class_<StandardThermoModelParamsNasa>(m, "StandardThermoModelParamsNasa")
        .def(py::init<>())
        .def_readwrite("dHf",         &StandardThermoModelParamsNasa::dHf)
        .def_readwrite("dH0",         &StandardThermoModelParamsNasa::dH0)
        .def_readwrite("H0",          &StandardThermoModelParamsNasa::H0)
        .def_readwrite("T0",          &StandardThermoModelParamsNasa::T0)
        .def_readwrite("polynomials", &StandardThermoModelParamsNasa::polynomials)
        ;

    m.def("StandardThermoModelNasa", StandardThermoModelNasa);
}
