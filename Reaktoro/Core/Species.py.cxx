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
#include <pybind11/stl.h>
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Core/FormationReaction.hpp>
#include <Reaktoro/Core/Species.hpp>
using namespace Reaktoro;

void exportSpecies(py::module& m)
{
    py::class_<Species>(m, "Species")
        .def(py::init<>())
        .def(py::init<String>())
        .def("clone", &Species::clone)
        .def("withName", &Species::withName)
        .def("withFormula", &Species::withFormula)
        .def("withSubstance", &Species::withSubstance)
        .def("withElements", &Species::withElements)
        .def("withCharge", &Species::withCharge)
        .def("withAggregateState", &Species::withAggregateState)
        .def("withFormationReaction", &Species::withFormationReaction)
        .def("withStandardGibbsEnergy", &Species::withStandardGibbsEnergy)
        .def("withStandardThermoModel", &Species::withStandardThermoModel)
        .def("withTags", &Species::withTags)
        .def("withAttachedData", &Species::withAttachedData)
        .def("name", &Species::name)
        .def("formula", &Species::formula)
        .def("substance", &Species::substance)
        .def("elements", &Species::elements, py::return_value_policy::reference_internal)
        .def("charge", &Species::charge)
        .def("aggregateState", &Species::aggregateState)
        .def("reaction", &Species::reaction, py::return_value_policy::reference_internal)
        .def("standardThermoModel", &Species::standardThermoModel, py::return_value_policy::reference_internal)
        .def("tags", &Species::tags, py::return_value_policy::reference_internal)
        .def("attachedData", &Species::attachedData)
        .def("molarMass", &Species::molarMass)
        .def("props", &Species::props)
        ;
}
