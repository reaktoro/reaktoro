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

#include "PyCompositionProblem.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktoro/Util/CompositionProblem.hpp>

namespace Reaktoro {

auto export_CompositionProblem() -> void
{
    py::class_<CompositionProblem>("CompositionProblem")
        .def(py::init<>())
        .def(py::init<const ChemicalSystem&>())
        .def("system", &CompositionProblem::system, py::return_internal_reference<>())
        .def("partition", &CompositionProblem::partition, py::return_internal_reference<>())
        .def("setPartition", &CompositionProblem::setPartition)
        .def("setTemperature", &CompositionProblem::setTemperature)
        .def("setPressure", &CompositionProblem::setPressure)
        .def("setAqueousFluid", &CompositionProblem::setAqueousComposition)
        .def("setGaseousFluid", &CompositionProblem::setGaseousComposition)
        .def("setSolid", &CompositionProblem::setSolidComposition)
        .def("setAqueousSaturation", &CompositionProblem::setAqueousSaturation)
        .def("setGaseousSaturation", &CompositionProblem::setGaseousSaturation)
        .def("setPorosity", &CompositionProblem::setPorosity)
        ;

    py::implicitly_convertible<CompositionProblem, EquilibriumInverseProblem>();
}

} // namespace Reaktoro
