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
#include <Reaktoro/Optimization/OptimumProblem.hpp>

namespace Reaktoro {

auto export_EquilibriumProblem() -> void
{
    auto setPartition1 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(const Partition&)>(&EquilibriumProblem::setPartition);
    auto setPartition2 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(std::string)>(&EquilibriumProblem::setPartition);

    auto setElementAmounts1 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(const Vector&)>(&EquilibriumProblem::setElementAmounts);
    auto setElementAmounts2 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(double)>(&EquilibriumProblem::setElementAmounts);

    auto setTemperature1 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(double)>(&EquilibriumProblem::setTemperature);
    auto setTemperature2 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(double, std::string)>(&EquilibriumProblem::setTemperature);

    auto setPressure1 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(double)>(&EquilibriumProblem::setPressure);
    auto setPressure2 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(double, std::string)>(&EquilibriumProblem::setPressure);

    auto add1 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(std::string, double, std::string)>(&EquilibriumProblem::add);
    auto add2 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(const ChemicalState&, double)>(&EquilibriumProblem::add);

    auto setSpeciesActivity1 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(std::string, double)>(&EquilibriumProblem::setSpeciesActivity);
    auto setSpeciesActivity2 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(std::string, double, std::string)>(&EquilibriumProblem::setSpeciesActivity);
    auto setSpeciesActivity3 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(std::string, double, std::string, std::string)>(&EquilibriumProblem::setSpeciesActivity);

    auto pH1 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(double)>(&EquilibriumProblem::pH);
    auto pH2 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(double, std::string)>(&EquilibriumProblem::pH);
    auto pH3 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(double, std::string, std::string)>(&EquilibriumProblem::pH);

    auto pe1 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(double)>(&EquilibriumProblem::pe);
    auto pe2 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(double, std::string)>(&EquilibriumProblem::pe);

    auto Eh1 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(double)>(&EquilibriumProblem::Eh);
    auto Eh2 = static_cast<EquilibriumProblem&(EquilibriumProblem::*)(double, std::string)>(&EquilibriumProblem::Eh);

    py::class_<EquilibriumProblem>("EquilibriumProblem", py::no_init)
        .def(py::init<const ChemicalSystem&>())
        .def(py::init<const EquilibriumProblem&>())
        .def("setPartition", setPartition1, py::return_internal_reference<>())
        .def("setPartition", setPartition2, py::return_internal_reference<>())
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
        .def("setSpeciesAmount", &EquilibriumProblem::setSpeciesAmount, py::return_internal_reference<>())
        .def("setSpeciesActivity", setSpeciesActivity1, py::return_internal_reference<>())
        .def("setSpeciesActivity", setSpeciesActivity2, py::return_internal_reference<>())
        .def("setSpeciesActivity", setSpeciesActivity3, py::return_internal_reference<>())
        .def("setPhaseAmount", &EquilibriumProblem::setPhaseAmount, py::return_internal_reference<>())
        .def("setPhaseVolume", &EquilibriumProblem::setPhaseVolume, py::return_internal_reference<>())
        .def("pH", pH1, py::return_internal_reference<>())
        .def("pH", pH2, py::return_internal_reference<>())
        .def("pH", pH3, py::return_internal_reference<>())
        .def("pe", pe1, py::return_internal_reference<>())
        .def("pe", pe2, py::return_internal_reference<>())
        .def("Eh", Eh1, py::return_internal_reference<>())
        .def("Eh", Eh2, py::return_internal_reference<>())
        .def("isInverseProblem", &EquilibriumProblem::isInverseProblem)
        .def("temperature", &EquilibriumProblem::temperature)
        .def("pressure", &EquilibriumProblem::pressure)
        .def("elementAmounts", &EquilibriumProblem::elementAmounts, py::return_value_policy<py::copy_const_reference>())
        .def("system", &EquilibriumProblem::system, py::return_value_policy<py::copy_const_reference>())
        .def("partition", &EquilibriumProblem::partition, py::return_value_policy<py::copy_const_reference>())
        ;
}

} // namespace Reaktoro
