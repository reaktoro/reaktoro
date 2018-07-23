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

#include "PyChemicalEditor.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>
#include <Reaktoro/Thermodynamics/Core/ChemicalEditor.hpp>
#include <Reaktoro/Thermodynamics/Core/Database.hpp>
#include <Reaktoro/Thermodynamics/Phases/AqueousPhase.hpp>
#include <Reaktoro/Thermodynamics/Phases/GaseousPhase.hpp>
#include <Reaktoro/Thermodynamics/Phases/MineralPhase.hpp>
#include <Reaktoro/Thermodynamics/Species/AqueousSpecies.hpp>
#include <Reaktoro/Thermodynamics/Species/GaseousSpecies.hpp>
#include <Reaktoro/Thermodynamics/Species/MineralSpecies.hpp>
#include <Reaktoro/Thermodynamics/Reactions/MineralReaction.hpp>

namespace Reaktoro {

auto export_ChemicalEditor() -> void
{
    using return_internal_ref = py::return_internal_reference<>;

    auto addPhase1 = static_cast<AqueousPhase&(ChemicalEditor::*)(const AqueousPhase&)>(&ChemicalEditor::addPhase);
    auto addPhase2 = static_cast<GaseousPhase&(ChemicalEditor::*)(const GaseousPhase&)>(&ChemicalEditor::addPhase);
    auto addPhase3 = static_cast<MineralPhase&(ChemicalEditor::*)(const MineralPhase&)>(&ChemicalEditor::addPhase);

    auto addAqueousPhase1 = static_cast<AqueousPhase&(ChemicalEditor::*)(std::vector<std::string>)>(&ChemicalEditor::addAqueousPhase);
    auto addAqueousPhase2 = static_cast<AqueousPhase&(ChemicalEditor::*)(std::string)>(&ChemicalEditor::addAqueousPhase);

    auto addGaseousPhase1 = static_cast<GaseousPhase&(ChemicalEditor::*)(std::vector<std::string>)>(&ChemicalEditor::addGaseousPhase);
    auto addGaseousPhase2 = static_cast<GaseousPhase&(ChemicalEditor::*)(std::string)>(&ChemicalEditor::addGaseousPhase);

    auto addMineralPhase1 = static_cast<MineralPhase&(ChemicalEditor::*)(std::vector<std::string>)>(&ChemicalEditor::addMineralPhase);
    auto addMineralPhase2 = static_cast<MineralPhase&(ChemicalEditor::*)(std::string)>(&ChemicalEditor::addMineralPhase);

    auto addMineralReaction1 = static_cast<MineralReaction&(ChemicalEditor::*)(const MineralReaction&)>(&ChemicalEditor::addMineralReaction);
    auto addMineralReaction2 = static_cast<MineralReaction&(ChemicalEditor::*)(std::string)>(&ChemicalEditor::addMineralReaction);

    auto aqueousPhase1 = static_cast<const AqueousPhase&(ChemicalEditor::*)() const>(&ChemicalEditor::aqueousPhase);
    auto aqueousPhase2 = static_cast<AqueousPhase&(ChemicalEditor::*)()>(&ChemicalEditor::aqueousPhase);

    auto gaseousPhase1 = static_cast<const GaseousPhase&(ChemicalEditor::*)() const>(&ChemicalEditor::gaseousPhase);
    auto gaseousPhase2 = static_cast<GaseousPhase&(ChemicalEditor::*)()>(&ChemicalEditor::gaseousPhase);

    auto mineralPhases1 = static_cast<const std::vector<MineralPhase>&(ChemicalEditor::*)() const>(&ChemicalEditor::mineralPhases);
    auto mineralPhases2 = static_cast<std::vector<MineralPhase>&(ChemicalEditor::*)()>(&ChemicalEditor::mineralPhases);

    py::class_<ChemicalEditor>("ChemicalEditor", py::no_init)
        .def(py::init<>())
        .def(py::init<const Database&>())
        .def("setTemperatures", &ChemicalEditor::setTemperatures)
        .def("setPressures", &ChemicalEditor::setPressures)
        .def("addPhase", addPhase1, return_internal_ref())
        .def("addPhase", addPhase2, return_internal_ref())
        .def("addPhase", addPhase3, return_internal_ref())
        .def("addReaction", &ChemicalEditor::addReaction, return_internal_ref())
        .def("addAqueousPhase", addAqueousPhase1, return_internal_ref())
        .def("addAqueousPhase", addAqueousPhase2, return_internal_ref())
        .def("addGaseousPhase", addGaseousPhase1, return_internal_ref())
        .def("addGaseousPhase", addGaseousPhase2, return_internal_ref())
        .def("addMineralPhase", addMineralPhase1, return_internal_ref())
        .def("addMineralPhase", addMineralPhase2, return_internal_ref())
        .def("addMineralReaction", addMineralReaction1, return_internal_ref())
        .def("addMineralReaction", addMineralReaction2, return_internal_ref())
        .def("aqueousPhase", aqueousPhase1, return_internal_ref())
        .def("aqueousPhase", aqueousPhase2, return_internal_ref())
        .def("gaseousPhase", gaseousPhase1, return_internal_ref())
        .def("gaseousPhase", gaseousPhase2, return_internal_ref())
        .def("mineralPhases", mineralPhases1, return_internal_ref())
        .def("mineralPhases", mineralPhases2, return_internal_ref())
        .def("createChemicalSystem", &ChemicalEditor::createChemicalSystem)
        .def("createReactionSystem", &ChemicalEditor::createReactionSystem)
        ;
}

} // namespace Reaktoro
