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
#include <Reaktoro/Water/WaterElectroProps.hpp>
using namespace Reaktoro;

void exportWaterElectroProps(py::module& m)
{
    py::class_<WaterElectroProps>(m, "WaterElectroProps")
        .def(py::init<>())
        .def_readwrite("epsilon", &WaterElectroProps::epsilon)
        .def_readwrite("epsilonT", &WaterElectroProps::epsilonT)
        .def_readwrite("epsilonP", &WaterElectroProps::epsilonP)
        .def_readwrite("epsilonTT", &WaterElectroProps::epsilonTT)
        .def_readwrite("epsilonTP", &WaterElectroProps::epsilonTP)
        .def_readwrite("epsilonPP", &WaterElectroProps::epsilonPP)
        .def_readwrite("bornZ", &WaterElectroProps::bornZ)
        .def_readwrite("bornY", &WaterElectroProps::bornY)
        .def_readwrite("bornQ", &WaterElectroProps::bornQ)
        .def_readwrite("bornN", &WaterElectroProps::bornN)
        .def_readwrite("bornU", &WaterElectroProps::bornU)
        .def_readwrite("bornX", &WaterElectroProps::bornX)
        ;
}
