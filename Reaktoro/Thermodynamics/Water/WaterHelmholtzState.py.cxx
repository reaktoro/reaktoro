// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Water/WaterHelmholtzState.hpp>
using namespace Reaktoro;

void exportWaterHelmholtzState(py::module& m)
{
    py::class_<WaterHelmholtzState>(m, "WaterHelmholtzState")
        .def(py::init<>())
        .def_readwrite("helmholtz", &WaterHelmholtzState::helmholtz)
        .def_readwrite("helmholtzT", &WaterHelmholtzState::helmholtzT)
        .def_readwrite("helmholtzD", &WaterHelmholtzState::helmholtzD)
        .def_readwrite("helmholtzTT", &WaterHelmholtzState::helmholtzTT)
        .def_readwrite("helmholtzTD", &WaterHelmholtzState::helmholtzTD)
        .def_readwrite("helmholtzDD", &WaterHelmholtzState::helmholtzDD)
        .def_readwrite("helmholtzTTT", &WaterHelmholtzState::helmholtzTTT)
        .def_readwrite("helmholtzTTD", &WaterHelmholtzState::helmholtzTTD)
        .def_readwrite("helmholtzTDD", &WaterHelmholtzState::helmholtzTDD)
        .def_readwrite("helmholtzDDD", &WaterHelmholtzState::helmholtzDDD)
        ;
}
