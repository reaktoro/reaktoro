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
#include <Reaktoro/Common/ChemicalScalar.hpp>
#include <Reaktoro/Common/ChemicalVector.hpp>
#include <Reaktoro/Common/ReactionEquation.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Reaction.hpp>

namespace Reaktoro {

void exportReaction(py::module& m)
{
    auto rate1 = static_cast<const ReactionRateFunction&(Reaction::*)() const>(&Reaction::rate);
    auto rate2 = static_cast<ChemicalScalar(Reaction::*)(const ChemicalProperties&) const>(&Reaction::rate);

    py::class_<Reaction>(m, "Reaction")
        .def(py::init<>())
        .def("setName", &Reaction::setName)
        .def("setEquilibriumConstant", &Reaction::setEquilibriumConstant)
        .def("setRate", &Reaction::setRate)
        .def("name", &Reaction::name)
        .def("equilibriumConstant", &Reaction::equilibriumConstant, py::return_value_policy::reference_internal)
        .def("equation", &Reaction::equation, py::return_value_policy::reference_internal)
        .def("system", &Reaction::system, py::return_value_policy::reference_internal)
        .def("species", &Reaction::species, py::return_value_policy::reference_internal)
        .def("indices", &Reaction::indices, py::return_value_policy::reference_internal)
        .def("stoichiometries", &Reaction::stoichiometries, py::return_value_policy::reference_internal)
        .def("stoichiometry", &Reaction::stoichiometry)
        .def("lnEquilibriumConstant", &Reaction::lnEquilibriumConstant)
        .def("lnReactionQuotient", &Reaction::lnReactionQuotient)
        .def("rate", rate1, py::return_value_policy::reference_internal)
        .def("rate", rate2)
        ;
}

} // namespace Reaktoro
