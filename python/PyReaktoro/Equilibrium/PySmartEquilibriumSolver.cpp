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
#include <Reaktoro/Equilibrium/SmartEquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumSolverClustering.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumSolverPriorityQueue.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumSolverNN.hpp>

namespace Reaktoro {

void exportSmartEquilibriumSolverClustering(py::module& m)
{
    auto solve1 = static_cast<SmartEquilibriumResult(SmartEquilibriumSolverClustering::*)(ChemicalState&, double, double, VectorConstRef)>(&SmartEquilibriumSolverClustering::solve);
    auto solve2 = static_cast<SmartEquilibriumResult(SmartEquilibriumSolverClustering::*)(ChemicalState&, const EquilibriumProblem&)>(&SmartEquilibriumSolverClustering::solve);

    py::class_<SmartEquilibriumSolverClustering>(m, "SmartEquilibriumSolverClustering")
        .def(py::init<const ChemicalSystem&>())
        .def(py::init<const Partition&>())
        .def("setOptions", &SmartEquilibriumSolverClustering::setOptions)
        .def("solve", solve1)
        .def("solve", solve2)
        .def("getProperties", &SmartEquilibriumSolverClustering::getProperties, py::return_value_policy::reference_internal)
        .def("getResult", &SmartEquilibriumSolverClustering::getResult, py::return_value_policy::reference_internal)
        .def("outputClusterInfo", &SmartEquilibriumSolverClustering::outputClusterInfo);
}

void exportSmartEquilibriumSolverPriorityQueue(py::module& m)
{
    auto solve1 = static_cast<SmartEquilibriumResult(SmartEquilibriumSolverPriorityQueue::*)(ChemicalState&, double, double, VectorConstRef)>(&SmartEquilibriumSolverPriorityQueue::solve);
    auto solve2 = static_cast<SmartEquilibriumResult(SmartEquilibriumSolverPriorityQueue::*)(ChemicalState&, const EquilibriumProblem&)>(&SmartEquilibriumSolverPriorityQueue::solve);

    py::class_<SmartEquilibriumSolverPriorityQueue>(m, "SmartEquilibriumSolverPriorityQueue")
        .def(py::init<const ChemicalSystem&>())
        .def(py::init<const Partition&>())
        .def("setOptions", &SmartEquilibriumSolverPriorityQueue::setOptions)
        .def("solve", solve1)
        .def("solve", solve2)
        .def("getProperties", &SmartEquilibriumSolverPriorityQueue::getProperties, py::return_value_policy::reference_internal)
        .def("getResult", &SmartEquilibriumSolverPriorityQueue::getResult, py::return_value_policy::reference_internal);
}

    void exportSmartEquilibriumSolverNN(py::module& m)
    {
        auto solve1 = static_cast<SmartEquilibriumResult(SmartEquilibriumSolverNN::*)(ChemicalState&, double, double, VectorConstRef)>(&SmartEquilibriumSolverNN::solve);
        auto solve2 = static_cast<SmartEquilibriumResult(SmartEquilibriumSolverNN::*)(ChemicalState&, const EquilibriumProblem&)>(&SmartEquilibriumSolverNN::solve);

        py::class_<SmartEquilibriumSolverNN>(m, "SmartEquilibriumSolverNN")
                .def(py::init<const ChemicalSystem&>())
                .def(py::init<const Partition&>())
                .def("setOptions", &SmartEquilibriumSolverNN::setOptions)
                .def("solve", solve1)
                .def("solve", solve2)
                .def("getProperties", &SmartEquilibriumSolverNN::getProperties, py::return_value_policy::reference_internal)
                .def("getResult", &SmartEquilibriumSolverNN::getResult, py::return_value_policy::reference_internal);
    }
} // namespace Reaktoro
