// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Reactions/MineralMechanism.hpp>

namespace Reaktoro {

void exportMineralMechanism(py::module& m)
{
    auto setCatalysts1 = static_cast<MineralMechanism&(MineralMechanism::*)(std::string)>(&MineralMechanism::setCatalysts);
    auto setCatalysts2 = static_cast<MineralMechanism&(MineralMechanism::*)(const MineralCatalyst&)>(&MineralMechanism::setCatalysts);
    auto setCatalysts3 = static_cast<MineralMechanism&(MineralMechanism::*)(const std::vector<MineralCatalyst>&)>(&MineralMechanism::setCatalysts);

    py::class_<MineralMechanism>(m, "MineralMechanism")
        .def(py::init<>())
        .def(py::init<std::string>())
        .def("setRateConstant", &MineralMechanism::setRateConstant, py::return_value_policy::reference_internal)
        .def("setActivationEnergy", &MineralMechanism::setActivationEnergy, py::return_value_policy::reference_internal)
        .def("setPowerP", &MineralMechanism::setPowerP, py::return_value_policy::reference_internal)
        .def("setPowerQ", &MineralMechanism::setPowerQ, py::return_value_policy::reference_internal)
        .def("setCatalysts", setCatalysts1, py::return_value_policy::reference_internal)
        .def("setCatalysts", setCatalysts2, py::return_value_policy::reference_internal)
        .def("setCatalysts", setCatalysts3, py::return_value_policy::reference_internal)
        .def_readwrite("kappa", &MineralMechanism::kappa)
        .def_readwrite("Ea", &MineralMechanism::Ea)
        .def_readwrite("p", &MineralMechanism::p)
        .def_readwrite("q", &MineralMechanism::q)
        .def_readwrite("catalysts", &MineralMechanism::catalysts)
        ;
}

} // namespace Reaktoro
