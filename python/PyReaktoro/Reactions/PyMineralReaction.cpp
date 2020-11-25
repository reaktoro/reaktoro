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

#include <PyReaktoro/PyReaktoro.hpp>

// Reaktoro includes
#include <Reaktoro/Common/ReactionEquation.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Reaction.hpp>
#include <Reaktoro/Thermodynamics/Reactions/MineralMechanism.hpp>
#include <Reaktoro/Thermodynamics/Reactions/MineralReaction.hpp>

namespace Reaktoro {

void exportMineralReaction(py::module& m)
{
    auto setEquation1 = static_cast<MineralReaction&(MineralReaction::*)(const ReactionEquation&)>(&MineralReaction::setEquation);
    auto setEquation2 = static_cast<MineralReaction&(MineralReaction::*)(std::string)>(&MineralReaction::setEquation);

    auto addMechanism1 = static_cast<MineralReaction&(MineralReaction::*)(const MineralMechanism&)>(&MineralReaction::addMechanism);
    auto addMechanism2 = static_cast<MineralReaction&(MineralReaction::*)(std::string)>(&MineralReaction::addMechanism);

    py::class_<MineralReaction>(m, "MineralReaction")
        .def(py::init<>())
        .def(py::init<std::string>())
        .def("setMineral", &MineralReaction::setMineral, py::return_value_policy::reference_internal)
        .def("setEquation", setEquation1, py::return_value_policy::reference_internal)
        .def("setEquation", setEquation2, py::return_value_policy::reference_internal)
        .def("setEquilibriumConstant", &MineralReaction::setEquilibriumConstant, py::return_value_policy::reference_internal)
        .def("setSpecificSurfaceArea", &MineralReaction::setSpecificSurfaceArea, py::return_value_policy::reference_internal)
        .def("setSurfaceArea", &MineralReaction::setSurfaceArea, py::return_value_policy::reference_internal)
        .def("addMechanism", addMechanism1, py::return_value_policy::reference_internal)
        .def("addMechanism", addMechanism2, py::return_value_policy::reference_internal)
        .def("setMechanisms", &MineralReaction::setMechanisms, py::return_value_policy::reference_internal)
        .def("mineral", &MineralReaction::mineral)
        .def("equation", &MineralReaction::equation, py::return_value_policy::reference_internal)
        .def("equilibriumConstant", &MineralReaction::equilibriumConstant, py::return_value_policy::reference_internal)
        .def("specificSurfaceArea", &MineralReaction::specificSurfaceArea)
        .def("volumetricSurfaceArea", &MineralReaction::volumetricSurfaceArea)
        .def("mechanisms", &MineralReaction::mechanisms, py::return_value_policy::reference_internal)
        ;

    m.def("createReaction", &createReaction);
}

} // namespace Reaktoro
