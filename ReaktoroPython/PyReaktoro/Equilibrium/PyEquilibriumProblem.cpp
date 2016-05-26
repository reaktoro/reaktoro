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

#include "PyEquilibriumProblem.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProblem.hpp>

namespace Reaktoro {

auto export_EquilibriumProblem() -> void
{
    auto setElementAmounts1 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(const Vector&)>(&EquilibriumProblem::setElementAmounts);
    auto setElementAmounts2 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(double)>(&EquilibriumProblem::setElementAmounts);

    auto setTemperature1 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(double)>(&EquilibriumProblem::setTemperature);
    auto setTemperature2 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(double, std::string)>(&EquilibriumProblem::setTemperature);

    auto setPressure1 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(double)>(&EquilibriumProblem::setPressure);
    auto setPressure2 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(double, std::string)>(&EquilibriumProblem::setPressure);

    auto add1 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(std::string, double, std::string)>(&EquilibriumProblem::add);
    auto add2 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(const ChemicalState&, double)>(&EquilibriumProblem::add);

    py::class_<EquilibriumProblem>("EquilibriumProblem", py::no_init)
        .def(py::init<const ChemicalSystem&>())
        .def(py::init<const EquilibriumProblem&>())
        .def("setPartition", &EquilibriumProblem::setPartition, py::return_internal_reference<>())
        .def("setTemperature", setTemperature1, py::return_internal_reference<>())
        .def("setTemperature", setTemperature2, py::return_internal_reference<>())
        .def("setPressure", setPressure1, py::return_internal_reference<>())
        .def("setPressure", setPressure2, py::return_internal_reference<>())
        .def("setElementAmounts", setElementAmounts1, py::return_internal_reference<>())
        .def("setElementAmounts", setElementAmounts2, py::return_internal_reference<>())
        .def("add", add1, py::return_internal_reference<>())
        .def("add", add2, py::return_internal_reference<>())
        .def("addCompound", &EquilibriumProblem::addCompound, py::return_internal_reference<>())
        .def("addSpecies", &EquilibriumProblem::addSpecies, py::return_internal_reference<>())
        .def("addState", &EquilibriumProblem::addState, py::return_internal_reference<>())
        .def("system", &EquilibriumProblem::system, py::return_internal_reference<>())
        .def("partition", &EquilibriumProblem::partition, py::return_internal_reference<>())
        .def("temperature", &EquilibriumProblem::temperature)
        .def("pressure", &EquilibriumProblem::pressure)
        .def("elementAmounts", &EquilibriumProblem::elementAmounts, py::return_internal_reference<>())
        ;
}

} // namespace Reaktoro
