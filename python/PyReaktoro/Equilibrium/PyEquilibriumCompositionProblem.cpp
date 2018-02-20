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
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumCompositionProblem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProblem.hpp>

namespace Reaktoro {

void exportEquilibriumCompositionProblem(py::module& m)
{
    // TODO remove this
//    py::class_<EquilibriumCompositionProblem>(m, "EquilibriumCompositionProblem")
//        .def(py::init<>())
//        .def(py::init<const ChemicalSystem&>())
//        .def("system", &EquilibriumCompositionProblem::system, py::return_value_policy::reference_internal)
//        .def("partition", &EquilibriumCompositionProblem::partition, py::return_value_policy::reference_internal)
//        .def("setPartition", &EquilibriumCompositionProblem::setPartition)
//        .def("setTemperature", &EquilibriumCompositionProblem::setTemperature)
//        .def("setPressure", &EquilibriumCompositionProblem::setPressure)
//        .def("setAqueousComposition", &EquilibriumCompositionProblem::setAqueousComposition)
//        .def("setGaseousComposition", &EquilibriumCompositionProblem::setGaseousComposition)
//        .def("setSolidComposition", &EquilibriumCompositionProblem::setSolidComposition)
//        .def("setAqueousSaturation", &EquilibriumCompositionProblem::setAqueousSaturation)
//        .def("setGaseousSaturation", &EquilibriumCompositionProblem::setGaseousSaturation)
//        .def("setPorosity", &EquilibriumCompositionProblem::setPorosity)
//        ;
//
//    py::implicitly_convertible<EquilibriumCompositionProblem, EquilibriumInverseProblem>();
}

} // namespace Reaktoro
