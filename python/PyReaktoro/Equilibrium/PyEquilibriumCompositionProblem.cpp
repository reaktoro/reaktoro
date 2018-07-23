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

#include "PyEquilibriumCompositionProblem.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumCompositionProblem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProblem.hpp>

namespace Reaktoro {

auto export_EquilibriumCompositionProblem() -> void
{
    py::class_<EquilibriumCompositionProblem>("EquilibriumCompositionProblem")
        .def(py::init<>())
        .def(py::init<const ChemicalSystem&>())
        .def("system", &EquilibriumCompositionProblem::system, py::return_internal_reference<>())
        .def("partition", &EquilibriumCompositionProblem::partition, py::return_internal_reference<>())
        .def("setPartition", &EquilibriumCompositionProblem::setPartition)
        .def("setTemperature", &EquilibriumCompositionProblem::setTemperature)
        .def("setPressure", &EquilibriumCompositionProblem::setPressure)
        .def("setAqueousComposition", &EquilibriumCompositionProblem::setAqueousComposition)
        .def("setGaseousComposition", &EquilibriumCompositionProblem::setGaseousComposition)
        .def("setSolidComposition", &EquilibriumCompositionProblem::setSolidComposition)
        .def("setAqueousSaturation", &EquilibriumCompositionProblem::setAqueousSaturation)
        .def("setGaseousSaturation", &EquilibriumCompositionProblem::setGaseousSaturation)
        .def("setPorosity", &EquilibriumCompositionProblem::setPorosity)
        ;

    py::implicitly_convertible<EquilibriumCompositionProblem, EquilibriumInverseProblem>();
}

} // namespace Reaktoro
