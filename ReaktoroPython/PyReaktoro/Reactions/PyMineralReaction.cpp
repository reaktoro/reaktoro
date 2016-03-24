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

#include "PyMineralReaction.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Common/ReactionEquation.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Reaction.hpp>
#include <Reaktoro/Thermodynamics/Reactions/MineralMechanism.hpp>
#include <Reaktoro/Thermodynamics/Reactions/MineralReaction.hpp>

namespace Reaktoro {

auto export_MineralReaction() -> void
{
    auto setEquation1 = static_cast<MineralReaction&(MineralReaction::*)(const ReactionEquation&)>(&MineralReaction::setEquation);
    auto setEquation2 = static_cast<MineralReaction&(MineralReaction::*)(std::string)>(&MineralReaction::setEquation);

    auto addMechanism1 = static_cast<MineralReaction&(MineralReaction::*)(const MineralMechanism&)>(&MineralReaction::addMechanism);
    auto addMechanism2 = static_cast<MineralReaction&(MineralReaction::*)(std::string)>(&MineralReaction::addMechanism);

    py::class_<MineralReaction>("MineralReaction")
        .def(py::init<>())
        .def(py::init<std::string>())
        .def("setMineral", &MineralReaction::setMineral, py::return_internal_reference<>())
        .def("setEquation", setEquation1, py::return_internal_reference<>())
        .def("setEquation", setEquation2, py::return_internal_reference<>())
        .def("setEquilibriumConstant", &MineralReaction::setEquilibriumConstant, py::return_internal_reference<>())
        .def("setSpecificSurfaceArea", &MineralReaction::setSpecificSurfaceArea, py::return_internal_reference<>())
        .def("addMechanism", addMechanism1, py::return_internal_reference<>())
        .def("addMechanism", addMechanism2, py::return_internal_reference<>())
        .def("setMechanisms", &MineralReaction::setMechanisms, py::return_internal_reference<>())
        .def("mineral", &MineralReaction::mineral)
        .def("equation", &MineralReaction::equation, py::return_internal_reference<>())
        .def("equilibriumConstant", &MineralReaction::equilibriumConstant, py::return_internal_reference<>())
        .def("specificSurfaceArea", &MineralReaction::specificSurfaceArea)
        .def("volumetricSurfaceArea", &MineralReaction::volumetricSurfaceArea)
        .def("mechanisms", &MineralReaction::mechanisms, py::return_internal_reference<>())
        ;
}

} // namespace Reaktoro
