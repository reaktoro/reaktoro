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
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProblem.hpp>

namespace Reaktoro {

void exportEquilibriumProblem(py::module& m)
{
    auto setElementAmounts1 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(VectorConstRef)>(&EquilibriumProblem::setElementAmounts);
    auto setElementAmounts2 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(double)>(&EquilibriumProblem::setElementAmounts);

    auto setElementAmount1 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(Index, double)>(&EquilibriumProblem::setElementAmount);
    auto setElementAmount2 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(std::string, double)>(&EquilibriumProblem::setElementAmount);

    auto setTemperature1 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(double)>(&EquilibriumProblem::setTemperature);
    auto setTemperature2 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(double, std::string)>(&EquilibriumProblem::setTemperature);

    auto setPressure1 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(double)>(&EquilibriumProblem::setPressure);
    auto setPressure2 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(double, std::string)>(&EquilibriumProblem::setPressure);

    auto add1 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(std::string, double, std::string)>(&EquilibriumProblem::add);
    auto add2 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(const ChemicalState&)>(&EquilibriumProblem::add);

    py::class_<EquilibriumProblem>(m, "EquilibriumProblem")
        .def(py::init<const ChemicalSystem&>())
        .def(py::init<const EquilibriumProblem&>())
        .def("setPartition", &EquilibriumProblem::setPartition, py::return_value_policy::reference_internal)
        .def("setTemperature", setTemperature1, py::return_value_policy::reference_internal)
        .def("setTemperature", setTemperature2, py::return_value_policy::reference_internal)
        .def("setPressure", setPressure1, py::return_value_policy::reference_internal)
        .def("setPressure", setPressure2, py::return_value_policy::reference_internal)
        .def("setElementAmounts", setElementAmounts1, py::return_value_policy::reference_internal)
        .def("setElementAmounts", setElementAmounts2, py::return_value_policy::reference_internal)
        .def("setElementAmount", setElementAmount1, py::return_value_policy::reference_internal)
        .def("setElementAmount", setElementAmount2, py::return_value_policy::reference_internal)
        .def("setElectricalCharge", &EquilibriumProblem::setElectricalCharge, py::return_value_policy::reference_internal)
        .def("add", add1, py::return_value_policy::reference_internal)
        .def("add", add2, py::return_value_policy::reference_internal)
        .def("addCompound", &EquilibriumProblem::addCompound, py::return_value_policy::reference_internal)
        .def("addSpecies", &EquilibriumProblem::addSpecies, py::return_value_policy::reference_internal)
        .def("addState", &EquilibriumProblem::addState, py::return_value_policy::reference_internal)
        .def("system", &EquilibriumProblem::system, py::return_value_policy::reference_internal)
        .def("partition", &EquilibriumProblem::partition, py::return_value_policy::reference_internal)
        .def("temperature", &EquilibriumProblem::temperature)
        .def("pressure", &EquilibriumProblem::pressure)
        .def("elementAmounts", &EquilibriumProblem::elementAmounts, py::return_value_policy::reference_internal)
        ;
}

} // namespace Reaktoro
