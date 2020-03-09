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

#include <PyReaktoro/PyReaktoro.hpp>

// Reaktoro includes
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProblem.hpp>

namespace Reaktoro {

void exportEquilibriumProblem(py::module& m)
{
    auto setElementAmount1 = static_cast<void(EquilibriumProblem::*)(Index, double)>(&EquilibriumProblem::setElementAmount);
    auto setElementAmount2 = static_cast<void(EquilibriumProblem::*)(std::string, double)>(&EquilibriumProblem::setElementAmount);

    auto setTemperature1 = static_cast<void(EquilibriumProblem::*)(double)>(&EquilibriumProblem::setTemperature);
    auto setTemperature2 = static_cast<void(EquilibriumProblem::*)(double, std::string)>(&EquilibriumProblem::setTemperature);

    auto setPressure1 = static_cast<void(EquilibriumProblem::*)(double)>(&EquilibriumProblem::setPressure);
    auto setPressure2 = static_cast<void(EquilibriumProblem::*)(double, std::string)>(&EquilibriumProblem::setPressure);

    auto add1 = static_cast<void(EquilibriumProblem::*)(std::string, double, std::string)>(&EquilibriumProblem::add);
    auto add2 = static_cast<void(EquilibriumProblem::*)(const ChemicalState&)>(&EquilibriumProblem::add);

    py::class_<EquilibriumProblem>(m, "EquilibriumProblem")
        .def(py::init<const ChemicalSystem&>())
        .def(py::init<const Partition&>())
        .def(py::init<const EquilibriumProblem&>())
        .def("setTemperature", setTemperature1)
        .def("setTemperature", setTemperature2)
        .def("setPressure", setPressure1)
        .def("setPressure", setPressure2)
        .def("setElementAmounts", &EquilibriumProblem::setElementAmounts)
        .def("setElementAmount", setElementAmount1)
        .def("setElementAmount", setElementAmount2)
        .def("setElectricalCharge", &EquilibriumProblem::setElectricalCharge)
        .def("add", add1)
        .def("add", add2)
        .def("addCompound", &EquilibriumProblem::addCompound)
        .def("addSpecies", &EquilibriumProblem::addSpecies)
        .def("addState", &EquilibriumProblem::addState)
        .def("system", &EquilibriumProblem::system, py::return_value_policy::reference_internal)
        .def("partition", &EquilibriumProblem::partition, py::return_value_policy::reference_internal)
        .def("temperature", &EquilibriumProblem::temperature)
        .def("pressure", &EquilibriumProblem::pressure)
        .def("elementAmounts", &EquilibriumProblem::elementAmounts, py::return_value_policy::reference_internal)

        // DEPRECATED METHODS: TO BE REMOVED
        .def("setPartition", &EquilibriumProblem::setPartition)
        ;
}

} // namespace Reaktoro
