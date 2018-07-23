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

#include "PyEquilibriumInverseProblem.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Equilibrium/EquilibriumInverseProblem.hpp>

namespace Reaktoro {

auto export_EquilibriumInverseProblem() -> void
{
    auto setTemperature1 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(double)>(&EquilibriumInverseProblem::setTemperature);
    auto setTemperature2 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(double, std::string)>(&EquilibriumInverseProblem::setTemperature);

    auto setPressure1 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(double)>(&EquilibriumInverseProblem::setPressure);
    auto setPressure2 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(double, std::string)>(&EquilibriumInverseProblem::setPressure);

    auto fixSpeciesAmount1 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(std::string, double, std::string)>(&EquilibriumInverseProblem::fixSpeciesAmount);
    auto fixSpeciesAmount2 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(std::string, double, std::string, std::string)>(&EquilibriumInverseProblem::fixSpeciesAmount);

    auto fixSpeciesMass1 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(std::string, double, std::string)>(&EquilibriumInverseProblem::fixSpeciesMass);
    auto fixSpeciesMass2 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(std::string, double, std::string, std::string)>(&EquilibriumInverseProblem::fixSpeciesMass);

    auto fixSpeciesActivity1 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(std::string, double)>(&EquilibriumInverseProblem::fixSpeciesActivity);
    auto fixSpeciesActivity2 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(std::string, double, std::string)>(&EquilibriumInverseProblem::fixSpeciesActivity);
    auto fixSpeciesActivity3 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(std::string, double, std::string, std::string)>(&EquilibriumInverseProblem::fixSpeciesActivity);

    auto fixSpeciesFugacity1 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(std::string, double, std::string)>(&EquilibriumInverseProblem::fixSpeciesFugacity);
    auto fixSpeciesFugacity2 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(std::string, double, std::string, std::string)>(&EquilibriumInverseProblem::fixSpeciesFugacity);

    auto pH1 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(double)>(&EquilibriumInverseProblem::pH);
    auto pH2 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(double, std::string)>(&EquilibriumInverseProblem::pH);
    auto pH3 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(double, std::string, std::string)>(&EquilibriumInverseProblem::pH);

    auto pE1 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(double)>(&EquilibriumInverseProblem::pE);
    auto pE2 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(double, std::string)>(&EquilibriumInverseProblem::pE);

    auto Eh1 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(double, std::string)>(&EquilibriumInverseProblem::Eh);
    auto Eh2 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(double, std::string, std::string)>(&EquilibriumInverseProblem::Eh);

    py::class_<EquilibriumInverseProblem>("EquilibriumInverseProblem", py::no_init)
        .def(py::init<const ChemicalSystem&>())
        .def(py::init<const EquilibriumInverseProblem&>())
        .def("setPartition", &EquilibriumInverseProblem::setPartition, py::return_internal_reference<>())
        .def("setTemperature", setTemperature1, py::return_internal_reference<>())
        .def("setTemperature", setTemperature2, py::return_internal_reference<>())
        .def("setPressure", setPressure1, py::return_internal_reference<>())
        .def("setPressure", setPressure2, py::return_internal_reference<>())
        .def("setElementInitialAmounts", &EquilibriumInverseProblem::setElementInitialAmounts, py::return_internal_reference<>())
        .def("add", &EquilibriumInverseProblem::add, py::return_internal_reference<>())
        .def("fixSpeciesAmount", fixSpeciesAmount1, py::return_internal_reference<>())
        .def("fixSpeciesAmount", fixSpeciesAmount2, py::return_internal_reference<>())
        .def("fixSpeciesMass", fixSpeciesMass1, py::return_internal_reference<>())
        .def("fixSpeciesMass", fixSpeciesMass2, py::return_internal_reference<>())
        .def("fixSpeciesActivity", fixSpeciesActivity1, py::return_internal_reference<>())
        .def("fixSpeciesActivity", fixSpeciesActivity2, py::return_internal_reference<>())
        .def("fixSpeciesActivity", fixSpeciesActivity3, py::return_internal_reference<>())
        .def("fixSpeciesFugacity", fixSpeciesFugacity1, py::return_internal_reference<>())
        .def("fixSpeciesFugacity", fixSpeciesFugacity2, py::return_internal_reference<>())
        .def("fixPhaseAmount", &EquilibriumInverseProblem::fixPhaseAmount, py::return_internal_reference<>())
        .def("fixPhaseMass", &EquilibriumInverseProblem::fixPhaseMass, py::return_internal_reference<>())
        .def("fixPhaseVolume", &EquilibriumInverseProblem::fixPhaseVolume, py::return_internal_reference<>())
		.def("fixPhaseSetVolume", &EquilibriumInverseProblem::fixPhaseSetVolume, py::return_internal_reference<>())
        .def("pH", pH1, py::return_internal_reference<>())
        .def("pH", pH2, py::return_internal_reference<>())
        .def("pH", pH3, py::return_internal_reference<>())
        .def("pE", pE1, py::return_internal_reference<>())
        .def("pE", pE2, py::return_internal_reference<>())
        .def("Eh", Eh1, py::return_internal_reference<>())
        .def("Eh", Eh2, py::return_internal_reference<>())
		.def("alkalinity", &EquilibriumInverseProblem::alkalinity, py::return_internal_reference<>())
        .def("system", &EquilibriumInverseProblem::system, py::return_internal_reference<>())
        .def("partition", &EquilibriumInverseProblem::partition, py::return_internal_reference<>())
        .def("temperature", &EquilibriumInverseProblem::temperature)
        .def("pressure", &EquilibriumInverseProblem::pressure)
        .def("numConstraints", &EquilibriumInverseProblem::numConstraints)
        .def("numTitrants", &EquilibriumInverseProblem::numTitrants)
        .def("formulaMatrixTitrants", &EquilibriumInverseProblem::formulaMatrixTitrants)
        .def("elementInitialAmounts", &EquilibriumInverseProblem::elementInitialAmounts)
        .def("titrantInitialAmounts", &EquilibriumInverseProblem::titrantInitialAmounts)
        ;
}

} // namespace Reaktoro
