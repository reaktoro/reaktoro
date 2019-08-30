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
#include <Reaktoro/Common/StringList.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>
#include <Reaktoro/Thermodynamics/Core/ChemicalEditor.hpp>
#include <Reaktoro/Thermodynamics/Core/Database.hpp>
#include <Reaktoro/Thermodynamics/Phases/AqueousPhase.hpp>
#include <Reaktoro/Thermodynamics/Phases/GaseousPhase.hpp>
#include <Reaktoro/Thermodynamics/Phases/LiquidPhase.hpp>
#include <Reaktoro/Thermodynamics/Phases/MineralPhase.hpp>
#include <Reaktoro/Thermodynamics/Reactions/MineralReaction.hpp>
#include <Reaktoro/Thermodynamics/Species/AqueousSpecies.hpp>
#include <Reaktoro/Thermodynamics/Species/GaseousSpecies.hpp>
#include <Reaktoro/Thermodynamics/Species/LiquidSpecies.hpp>
#include <Reaktoro/Thermodynamics/Species/MineralSpecies.hpp>

namespace Reaktoro {

void exportChemicalEditor(py::module& m)
{
    auto addPhase1 = static_cast<AqueousPhase&(ChemicalEditor::*)(const AqueousPhase&)>(&ChemicalEditor::addPhase);
    auto addPhase2 = static_cast<GaseousPhase&(ChemicalEditor::*)(const GaseousPhase&)>(&ChemicalEditor::addPhase);
    auto addPhase3 = static_cast<LiquidPhase&(ChemicalEditor::*)(const LiquidPhase&)>(&ChemicalEditor::addPhase);
    auto addPhase4 = static_cast<MineralPhase&(ChemicalEditor::*)(const MineralPhase&)>(&ChemicalEditor::addPhase);

    auto addMineralReaction1 = static_cast<MineralReaction&(ChemicalEditor::*)(const MineralReaction&)>(&ChemicalEditor::addMineralReaction);
    auto addMineralReaction2 = static_cast<MineralReaction&(ChemicalEditor::*)(std::string)>(&ChemicalEditor::addMineralReaction);

    auto aqueousPhase1 = static_cast<const AqueousPhase&(ChemicalEditor::*)() const>(&ChemicalEditor::aqueousPhase);
    auto aqueousPhase2 = static_cast<AqueousPhase&(ChemicalEditor::*)()>(&ChemicalEditor::aqueousPhase);

    auto gaseousPhase1 = static_cast<const GaseousPhase&(ChemicalEditor::*)() const>(&ChemicalEditor::gaseousPhase);
    auto gaseousPhase2 = static_cast<GaseousPhase&(ChemicalEditor::*)()>(&ChemicalEditor::gaseousPhase);

    auto liquidPhase1 = static_cast<const LiquidPhase&(ChemicalEditor::*)() const>(&ChemicalEditor::liquidPhase);
    auto liquidPhase2 = static_cast<LiquidPhase&(ChemicalEditor::*)()>(&ChemicalEditor::liquidPhase);

    auto mineralPhases1 = static_cast<const std::vector<MineralPhase>&(ChemicalEditor::*)() const>(&ChemicalEditor::mineralPhases);
    auto mineralPhases2 = static_cast<std::vector<MineralPhase>&(ChemicalEditor::*)()>(&ChemicalEditor::mineralPhases);

    py::class_<ChemicalEditor>(m, "ChemicalEditor")
        .def(py::init<>())
        .def(py::init<const Database&>())
        .def("setTemperatures", &ChemicalEditor::setTemperatures)
        .def("setPressures", &ChemicalEditor::setPressures)
        .def("addPhase", addPhase1, py::return_value_policy::reference_internal)
        .def("addPhase", addPhase2, py::return_value_policy::reference_internal)
        .def("addPhase", addPhase3, py::return_value_policy::reference_internal)
        .def("addPhase", addPhase4, py::return_value_policy::reference_internal)
        .def("addReaction", &ChemicalEditor::addReaction, py::return_value_policy::reference_internal)
        
        .def("addAqueousPhase", &ChemicalEditor::addAqueousPhase, py::return_value_policy::reference_internal)
        .def("addAqueousPhaseWithElements", &ChemicalEditor::addAqueousPhaseWithElements, py::return_value_policy::reference_internal)
        .def("addAqueousPhaseWithElementsOf", &ChemicalEditor::addAqueousPhaseWithElementsOf, py::return_value_policy::reference_internal)

        .def("addGaseousPhase", &ChemicalEditor::addGaseousPhase, py::return_value_policy::reference_internal)
        .def("addGaseousPhaseWithElements", &ChemicalEditor::addGaseousPhaseWithElements, py::return_value_policy::reference_internal)
        .def("addGaseousPhaseWithElementsOf", &ChemicalEditor::addGaseousPhaseWithElementsOf, py::return_value_policy::reference_internal)

        .def("addLiquidPhase", &ChemicalEditor::addLiquidPhase, py::return_value_policy::reference_internal)
        .def("addLiquidPhaseWithElements", &ChemicalEditor::addLiquidPhaseWithElements, py::return_value_policy::reference_internal)
        .def("addLiquidPhaseWithElementsOf", &ChemicalEditor::addLiquidPhaseWithElementsOf, py::return_value_policy::reference_internal)

        .def("addMineralPhase", &ChemicalEditor::addMineralPhase, py::return_value_policy::reference_internal)
        .def("addMineralPhaseWithElements", &ChemicalEditor::addMineralPhaseWithElements, py::return_value_policy::reference_internal)
        .def("addMineralPhaseWithElementsOf", &ChemicalEditor::addMineralPhaseWithElementsOf, py::return_value_policy::reference_internal)

        .def("addMineralReaction", addMineralReaction1, py::return_value_policy::reference_internal)
        .def("addMineralReaction", addMineralReaction2, py::return_value_policy::reference_internal)
        .def("aqueousPhase", aqueousPhase1, py::return_value_policy::reference_internal)
        .def("aqueousPhase", aqueousPhase2, py::return_value_policy::reference_internal)
        .def("gaseousPhase", gaseousPhase1, py::return_value_policy::reference_internal)
        .def("gaseousPhase", gaseousPhase2, py::return_value_policy::reference_internal)
        .def("liquidPhase", liquidPhase1, py::return_value_policy::reference_internal)
        .def("liquidPhase", liquidPhase2, py::return_value_policy::reference_internal)
        .def("mineralPhases", mineralPhases1, py::return_value_policy::reference_internal)
        .def("mineralPhases", mineralPhases2, py::return_value_policy::reference_internal)
        .def("createChemicalSystem", &ChemicalEditor::createChemicalSystem)
        .def("createReactionSystem", &ChemicalEditor::createReactionSystem)
        ;
}

} // namespace Reaktoro
