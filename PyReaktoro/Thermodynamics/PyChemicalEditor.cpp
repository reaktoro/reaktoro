// Reaktoro is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
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

#include "PyChemicalEditor.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>
#include <Reaktoro/Thermodynamics/Core/ChemicalEditor.hpp>
#include <Reaktoro/Thermodynamics/Core/Database.hpp>
#include <Reaktoro/Thermodynamics/Phase/AqueousPhase.hpp>
#include <Reaktoro/Thermodynamics/Phase/GaseousPhase.hpp>
#include <Reaktoro/Thermodynamics/Phase/MineralPhase.hpp>
#include <Reaktoro/Thermodynamics/Species/AqueousSpecies.hpp>
#include <Reaktoro/Thermodynamics/Species/GaseousSpecies.hpp>
#include <Reaktoro/Thermodynamics/Species/MineralSpecies.hpp>
#include <Reaktoro/Thermodynamics/Reactions/MineralReaction.hpp>

namespace Reaktoro {

auto export_ChemicalEditor() -> void
{
    using return_const_ref = py::return_value_policy<py::copy_const_reference>;

    auto addPhase1 = static_cast<AqueousPhase&(ChemicalEditor::*)(const AqueousPhase&)>(&ChemicalEditor::addPhase);
    auto addPhase2 = static_cast<GaseousPhase&(ChemicalEditor::*)(const GaseousPhase&)>(&ChemicalEditor::addPhase);
    auto addPhase3 = static_cast<MineralPhase&(ChemicalEditor::*)(const MineralPhase&)>(&ChemicalEditor::addPhase);

    auto addAqueousPhase1 = static_cast<AqueousPhase&(ChemicalEditor::*)(const std::vector<std::string>&)>(&ChemicalEditor::addAqueousPhase);
    auto addAqueousPhase2 = static_cast<AqueousPhase&(ChemicalEditor::*)(const std::string&)>(&ChemicalEditor::addAqueousPhase);
    
    auto addGaseousPhase1 = static_cast<GaseousPhase&(ChemicalEditor::*)(const std::vector<std::string>&)>(&ChemicalEditor::addGaseousPhase);
    auto addGaseousPhase2 = static_cast<GaseousPhase&(ChemicalEditor::*)(const std::string&)>(&ChemicalEditor::addGaseousPhase);

    auto addMineralPhase1 = static_cast<MineralPhase&(ChemicalEditor::*)(const std::vector<std::string>&)>(&ChemicalEditor::addMineralPhase);
    auto addMineralPhase2 = static_cast<MineralPhase&(ChemicalEditor::*)(const std::string&)>(&ChemicalEditor::addMineralPhase);

    auto aqueousPhase1 = static_cast<const AqueousPhase&(ChemicalEditor::*)() const>(&ChemicalEditor::aqueousPhase);
    auto aqueousPhase2 = static_cast<AqueousPhase&(ChemicalEditor::*)()>(&ChemicalEditor::aqueousPhase);

    auto gaseousPhase1 = static_cast<const GaseousPhase&(ChemicalEditor::*)() const>(&ChemicalEditor::gaseousPhase);
    auto gaseousPhase2 = static_cast<GaseousPhase&(ChemicalEditor::*)()>(&ChemicalEditor::gaseousPhase);

    auto mineralPhases1 = static_cast<const std::vector<MineralPhase>&(ChemicalEditor::*)() const>(&ChemicalEditor::mineralPhases);
    auto mineralPhases2 = static_cast<std::vector<MineralPhase>&(ChemicalEditor::*)()>(&ChemicalEditor::mineralPhases);

    py::class_<ChemicalEditor>("ChemicalEditor", py::no_init)
        .def(py::init<const Database&>())
        .def("setTemperatures", &ChemicalEditor::setTemperatures)
        .def("setPressures", &ChemicalEditor::setPressures)
        .def("addPhase", addPhase1)
        .def("addPhase", addPhase2)
        .def("addPhase", addPhase3)
        .def("addReaction", &ChemicalEditor::addReaction)
        .def("addAqueousPhase", addAqueousPhase1)
        .def("addAqueousPhase", addAqueousPhase2)
        .def("addGaseousPhase", addGaseousPhase1)
        .def("addGaseousPhase", addGaseousPhase2)
        .def("addMineralPhase", addMineralPhase1)
        .def("addMineralPhase", addMineralPhase2)
        .def("addMineralReaction", &ChemicalEditor::addMineralReaction)
        .def("aqueousPhase", aqueousPhase1)
        .def("aqueousPhase", aqueousPhase2)
        .def("gaseousPhase", gaseousPhase1)
        .def("gaseousPhase", gaseousPhase2)
        .def("mineralPhases", mineralPhases1)
        .def("mineralPhases", mineralPhases2)
        .def("createChemicalSystem", &ChemicalEditor::createChemicalSystem)
        .def("createReactionSystem", &ChemicalEditor::createReactionSystem)
        ;
}

} // namespace Reaktoro
