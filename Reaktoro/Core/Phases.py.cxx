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
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Core/Phases.hpp>
using namespace Reaktoro;

void exportPhases(py::module& m)
{
    py::class_<Speciate>(m, "Speciate")
        .def(py::init<>())
        .def_readwrite("symbols", &Speciate::symbols)
        ;

    m.def("speciate", speciate);

    py::class_<Phases>(m, "Phases")
        .def(py::init<const Database&>())
        .def("add", py::overload_cast<const GenericPhase&>(&Phases::add))
        .def("add", py::overload_cast<const GenericPhasesGenerator&>(&Phases::add))
        .def("database", &Phases::database)
        ;

    py::class_<GenericPhase>(m, "GenericPhase")
        .def(py::init<>())
        .def(py::init<const StringList&>())
        .def(py::init<const Speciate&>())
        .def("setName", &GenericPhase::setName, py::return_value_policy::reference_internal)
        .def("setStateOfMatter", &GenericPhase::setStateOfMatter, py::return_value_policy::reference_internal)
        .def("setAggregateState", &GenericPhase::setAggregateState, py::return_value_policy::reference_internal)
        .def("setActivityModel", &GenericPhase::setActivityModel, py::return_value_policy::reference_internal)
        .def("setIdealActivityModel", &GenericPhase::setIdealActivityModel, py::return_value_policy::reference_internal)
        .def("named", &GenericPhase::named, py::return_value_policy::reference_internal)
        .def("set", py::overload_cast<StateOfMatter>(&GenericPhase::set), py::return_value_policy::reference_internal)
        .def("set", py::overload_cast<AggregateState>(&GenericPhase::set), py::return_value_policy::reference_internal)
        .def("set", py::overload_cast<const ActivityModel&>(&GenericPhase::set), py::return_value_policy::reference_internal)
        .def("name", &GenericPhase::name)
        .def("stateOfMatter", &GenericPhase::stateOfMatter)
        .def("aggregateState", &GenericPhase::aggregateState)
        .def("species", &GenericPhase::species, py::return_value_policy::reference_internal)
        .def("elements", &GenericPhase::elements, py::return_value_policy::reference_internal)
        .def("activityModel", &GenericPhase::activityModel, py::return_value_policy::reference_internal)
        .def("idealActivityModel", &GenericPhase::idealActivityModel, py::return_value_policy::reference_internal)
        .def("convert", &GenericPhase::convert)
        ;

    py::class_<GenericPhasesGenerator>(m, "GenericPhasesGenerator")
        .def(py::init<>())
        .def(py::init<const StringList&>())
        .def(py::init<const Speciate&>())
        .def("setStateOfMatter", &GenericPhasesGenerator::setStateOfMatter, py::return_value_policy::reference_internal)
        .def("setAggregateState", &GenericPhasesGenerator::setAggregateState, py::return_value_policy::reference_internal)
        .def("setActivityModel", &GenericPhasesGenerator::setActivityModel, py::return_value_policy::reference_internal)
        .def("setIdealActivityModel", &GenericPhasesGenerator::setIdealActivityModel, py::return_value_policy::reference_internal)
        .def("set", py::overload_cast<StateOfMatter>(&GenericPhasesGenerator::set), py::return_value_policy::reference_internal)
        .def("set", py::overload_cast<AggregateState>(&GenericPhasesGenerator::set), py::return_value_policy::reference_internal)
        .def("set", py::overload_cast<const ActivityModel&>(&GenericPhasesGenerator::set), py::return_value_policy::reference_internal)
        .def("stateOfMatter", &GenericPhasesGenerator::stateOfMatter)
        .def("aggregateState", &GenericPhasesGenerator::aggregateState)
        .def("species", &GenericPhasesGenerator::species, py::return_value_policy::reference_internal)
        .def("elements", &GenericPhasesGenerator::elements, py::return_value_policy::reference_internal)
        .def("activityModel", &GenericPhasesGenerator::activityModel, py::return_value_policy::reference_internal)
        .def("idealActivityModel", &GenericPhasesGenerator::idealActivityModel, py::return_value_policy::reference_internal)
        .def("convert", &GenericPhasesGenerator::convert)
        ;

    py::class_<AqueousPhase, GenericPhase>(m, "AqueousPhase")
        .def(py::init<>())
        .def(py::init<const StringList&>())
        .def(py::init<Speciate>())
        ;

    py::class_<GaseousPhase, GenericPhase>(m, "GaseousPhase")
        .def(py::init<>())
        .def(py::init<const StringList&>())
        .def(py::init<const Speciate&>())
        ;

    py::class_<LiquidPhase, GenericPhase>(m, "LiquidPhase")
        .def(py::init<>())
        .def(py::init<const StringList&>())
        .def(py::init<const Speciate&>())
        ;

    py::class_<SolidPhase, GenericPhase>(m, "SolidPhase")
        .def(py::init<const StringList&>())
        ;

    py::class_<MineralPhase, GenericPhase>(m, "MineralPhase")
        .def(py::init<String>())
        ;

    py::class_<MineralPhases, GenericPhasesGenerator>(m, "MineralPhases")
        .def(py::init<>())
        .def(py::init<const StringList&>())
        .def(py::init<const Speciate&>())
        ;
}
