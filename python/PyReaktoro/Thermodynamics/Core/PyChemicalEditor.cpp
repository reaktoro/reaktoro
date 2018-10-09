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

// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

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

void exportChemicalEditor(py::module& m)
{
    auto addPhase1 = static_cast<AqueousPhase&(ChemicalEditor::*)(const AqueousPhase&)>(&ChemicalEditor::addPhase);
    auto addPhase2 = static_cast<GaseousPhase&(ChemicalEditor::*)(const GaseousPhase&)>(&ChemicalEditor::addPhase);
    auto addPhase3 = static_cast<MineralPhase&(ChemicalEditor::*)(const MineralPhase&)>(&ChemicalEditor::addPhase);

    //Need to be this order so that pybind do not see a python string as C++ std::string not a C++ std::vector<std::string> 
    //with each leter as an elemento of that vetor 
    auto addAqueousPhase1 = static_cast<AqueousPhase&(ChemicalEditor::*)(std::string)>(&ChemicalEditor::addAqueousPhase);
    auto addAqueousPhase2 = static_cast<AqueousPhase&(ChemicalEditor::*)(std::vector<std::string>)>(&ChemicalEditor::addAqueousPhase);

    //Need to be this order so that pybind do not see a python string as C++ std::string not a C++ std::vector<std::string> 
    //with each leter as an elemento of that vetor 
    auto addGaseousPhase1 = static_cast<GaseousPhase&(ChemicalEditor::*)(std::string)>(&ChemicalEditor::addGaseousPhase);
    auto addGaseousPhase2 = static_cast<GaseousPhase&(ChemicalEditor::*)(std::vector<std::string>)>(&ChemicalEditor::addGaseousPhase);

    //Need to be this order so that pybind do not see a python string as C++ std::string not a C++ std::vector<std::string> 
    //with each leter as an elemento of that vetor 
    auto addMineralPhase1 = static_cast<MineralPhase&(ChemicalEditor::*)(std::string)>(&ChemicalEditor::addMineralPhase);
    auto addMineralPhase2 = static_cast<MineralPhase&(ChemicalEditor::*)(std::vector<std::string>)>(&ChemicalEditor::addMineralPhase);

    auto addMineralReaction1 = static_cast<MineralReaction&(ChemicalEditor::*)(const MineralReaction&)>(&ChemicalEditor::addMineralReaction);
    auto addMineralReaction2 = static_cast<MineralReaction&(ChemicalEditor::*)(std::string)>(&ChemicalEditor::addMineralReaction);

    auto aqueousPhase1 = static_cast<const AqueousPhase&(ChemicalEditor::*)() const>(&ChemicalEditor::aqueousPhase);
    auto aqueousPhase2 = static_cast<AqueousPhase&(ChemicalEditor::*)()>(&ChemicalEditor::aqueousPhase);

    auto gaseousPhase1 = static_cast<const GaseousPhase&(ChemicalEditor::*)() const>(&ChemicalEditor::gaseousPhase);
    auto gaseousPhase2 = static_cast<GaseousPhase&(ChemicalEditor::*)()>(&ChemicalEditor::gaseousPhase);

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
        .def("addReaction", &ChemicalEditor::addReaction, py::return_value_policy::reference_internal)
        .def("addAqueousPhase", addAqueousPhase1, py::return_value_policy::reference_internal)
        .def("addAqueousPhase", addAqueousPhase2, py::return_value_policy::reference_internal)
        .def("addGaseousPhase", addGaseousPhase1, py::return_value_policy::reference_internal)
        .def("addGaseousPhase", addGaseousPhase2, py::return_value_policy::reference_internal)
        .def("addMineralPhase", addMineralPhase1, py::return_value_policy::reference_internal)
        .def("addMineralPhase", addMineralPhase2, py::return_value_policy::reference_internal)
        .def("addMineralReaction", addMineralReaction1, py::return_value_policy::reference_internal)
        .def("addMineralReaction", addMineralReaction2, py::return_value_policy::reference_internal)
        .def("aqueousPhase", aqueousPhase1, py::return_value_policy::reference_internal)
        .def("aqueousPhase", aqueousPhase2, py::return_value_policy::reference_internal)
        .def("gaseousPhase", gaseousPhase1, py::return_value_policy::reference_internal)
        .def("gaseousPhase", gaseousPhase2, py::return_value_policy::reference_internal)
        .def("mineralPhases", mineralPhases1, py::return_value_policy::reference_internal)
        .def("mineralPhases", mineralPhases2, py::return_value_policy::reference_internal)
        .def("createChemicalSystem", &ChemicalEditor::createChemicalSystem)
        .def("createReactionSystem", &ChemicalEditor::createReactionSystem)
        ;
}

} // namespace Reaktoro
