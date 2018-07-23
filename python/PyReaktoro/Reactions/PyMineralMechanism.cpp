// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

#include "PyMineralMechanism.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Reactions/MineralMechanism.hpp>

// PyReaktoro includes
#include <PyReaktoro/Common/PyConverters.hpp>

namespace Reaktoro {

auto export_MineralMechanism() -> void
{
    using return_internal_ref = py::return_internal_reference<>;

    auto setCatalysts1 = static_cast<MineralMechanism&(MineralMechanism::*)(std::string)>(&MineralMechanism::setCatalysts);
    auto setCatalysts2 = static_cast<MineralMechanism&(MineralMechanism::*)(const MineralCatalyst&)>(&MineralMechanism::setCatalysts);
    auto setCatalysts3 = static_cast<MineralMechanism&(MineralMechanism::*)(const std::vector<MineralCatalyst>&)>(&MineralMechanism::setCatalysts);

    py::class_<MineralMechanism>("MineralMechanism")
        .def(py::init<>())
        .def(py::init<std::string>())
        .def("setRateConstant", &MineralMechanism::setRateConstant, return_internal_ref())
        .def("setActivationEnergy", &MineralMechanism::setActivationEnergy, return_internal_ref())
        .def("setPowerP", &MineralMechanism::setPowerP, return_internal_ref())
        .def("setPowerQ", &MineralMechanism::setPowerQ, return_internal_ref())
        .def("setCatalysts", setCatalysts1, return_internal_ref())
        .def("setCatalysts", setCatalysts2, return_internal_ref())
        .def("setCatalysts", setCatalysts3, return_internal_ref())
        .def_readwrite("kappa", &MineralMechanism::kappa)
        .def_readwrite("Ea", &MineralMechanism::Ea)
        .def_readwrite("p", &MineralMechanism::p)
        .def_readwrite("q", &MineralMechanism::q)
        .def_readwrite("catalysts", &MineralMechanism::catalysts)
        ;

    export_std_vector<MineralMechanism>("MineralMechanismVector");
}

} // namespace Reaktoro
