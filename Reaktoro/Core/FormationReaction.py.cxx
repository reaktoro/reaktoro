// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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
#include <pybind11/functional.h>
#include <pybind11/stl.h>
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Core/FormationReaction.hpp>
#include <Reaktoro/Core/Species.hpp>
using namespace Reaktoro;

void exportFormationReaction(py::module& m)
{
    py::class_<FormationReaction>(m, "FormationReaction")
        .def(py::init<>())
        .def("clone", &FormationReaction::clone)
        .def("withReactants", &FormationReaction::withReactants)
        .def("withEquilibriumConstant", &FormationReaction::withEquilibriumConstant)
        .def("withProductStandardVolume", &FormationReaction::withProductStandardVolume)
        .def("withProductStandardVolumeModel", &FormationReaction::withProductStandardVolumeModel)
        .def("withReactionThermoModel", &FormationReaction::withReactionThermoModel)
        .def("reactants", &FormationReaction::reactants)
        .def("stoichiometry", &FormationReaction::stoichiometry)
        .def("productStandardVolumeModel", &FormationReaction::productStandardVolumeModel)
        .def("reactionThermoModel", &FormationReaction::reactionThermoModel)
        .def("createStandardThermoModel", &FormationReaction::createStandardThermoModel)
        ;
}
