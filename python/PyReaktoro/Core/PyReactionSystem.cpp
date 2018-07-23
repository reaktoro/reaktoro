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

#include "PyReactionSystem.hpp"

// Boost includes
#include <boost/python.hpp>
#include <boost/smart_ptr.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Common/ChemicalVector.hpp>
#include <Reaktoro/Common/ThermoVector.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Reaction.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>
#include <Reaktoro/Thermodynamics/Core/ChemicalEditor.hpp>

namespace Reaktoro {
namespace {

auto createReactionSystemFromChemicalEditor(const ChemicalEditor& editor) -> boost::shared_ptr<ReactionSystem>
{
    return boost::make_shared<ReactionSystem>(editor);
}

} // namespace

auto export_ReactionSystem() -> void
{
    auto reaction1 = static_cast<const Reaction&(ReactionSystem::*)(Index) const>(&ReactionSystem::reaction);
    auto reaction2 = static_cast<const Reaction&(ReactionSystem::*)(std::string) const>(&ReactionSystem::reaction);

    py::class_<ReactionSystem>("ReactionSystem")
        .def(py::init<>())
        .def(py::init<const ChemicalSystem&, const std::vector<Reaction>&>())
        .def("__init__", py::make_constructor(createReactionSystemFromChemicalEditor))
        .def("numReactions", &ReactionSystem::numReactions)
        .def("indexReaction", &ReactionSystem::indexReaction)
        .def("reactions", &ReactionSystem::reactions, py::return_internal_reference<>())
        .def("reaction", reaction1, py::return_internal_reference<>())
        .def("reaction", reaction2, py::return_internal_reference<>())
        .def("stoichiometricMatrix", &ReactionSystem::stoichiometricMatrix, py::return_internal_reference<>())
        .def("system", &ReactionSystem::system, py::return_internal_reference<>())
        .def("lnEquilibriumConstants", &ReactionSystem::lnEquilibriumConstants)
        .def("lnReactionQuotients", &ReactionSystem::lnReactionQuotients)
        .def("rates", &ReactionSystem::rates)
        ;
}

} // namespace Reaktoro
