// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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
        .def("withProduct", &FormationReaction::withProduct)
        .def("withReactants", &FormationReaction::withReactants)
        .def("withEquilibriumConstant", &FormationReaction::withEquilibriumConstant)
        .def("withEquilibriumConstantFn", &FormationReaction::withEquilibriumConstantFn)
        .def("withEnthalpyChange", &FormationReaction::withEnthalpyChange)
        .def("withEnthalpyChangeFn", &FormationReaction::withEnthalpyChangeFn)
        .def("product", &FormationReaction::product)
        .def("reactants", &FormationReaction::reactants)
        .def("equilibriumConstantFn", &FormationReaction::equilibriumConstantFn)
        .def("enthalpyChangeFn", &FormationReaction::enthalpyChangeFn)
        .def("standardGibbsEnergyFn", &FormationReaction::standardGibbsEnergyFn)
        .def("standardEnthalpyFn", &FormationReaction::standardEnthalpyFn)
        .def("stoichiometry", &FormationReaction::stoichiometry)
        ;
}
