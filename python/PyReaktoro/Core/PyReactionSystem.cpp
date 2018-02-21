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
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Common/ChemicalVector.hpp>
#include <Reaktoro/Common/ThermoVector.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Reaction.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>
#include <Reaktoro/Thermodynamics/Core/ChemicalEditor.hpp>

namespace Reaktoro {

void exportReactionSystem(py::module& m)
{
    auto reaction1 = static_cast<const Reaction&(ReactionSystem::*)(Index) const>(&ReactionSystem::reaction);
    auto reaction2 = static_cast<const Reaction&(ReactionSystem::*)(std::string) const>(&ReactionSystem::reaction);

    py::class_<ReactionSystem>(m, "ReactionSystem")
        .def(py::init<>())
        .def(py::init<const ChemicalSystem&, const std::vector<Reaction>&>())
        .def(py::init([](const ChemicalEditor& editor) { return std::make_unique<ReactionSystem>(editor); }))
        .def("numReactions", &ReactionSystem::numReactions)
        .def("indexReaction", &ReactionSystem::indexReaction)
        .def("reactions", &ReactionSystem::reactions, py::return_value_policy::reference_internal)
        .def("reaction", reaction1, py::return_value_policy::reference_internal)
        .def("reaction", reaction2, py::return_value_policy::reference_internal)
        .def("stoichiometricMatrix", &ReactionSystem::stoichiometricMatrix, py::return_value_policy::reference_internal)
        .def("system", &ReactionSystem::system, py::return_value_policy::reference_internal)
        .def("lnEquilibriumConstants", &ReactionSystem::lnEquilibriumConstants)
        .def("lnReactionQuotients", &ReactionSystem::lnReactionQuotients)
        .def("rates", &ReactionSystem::rates)
        ;
}

} // namespace Reaktoro
