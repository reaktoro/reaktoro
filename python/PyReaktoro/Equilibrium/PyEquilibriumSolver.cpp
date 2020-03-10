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
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>

namespace Reaktoro {

void exportEquilibriumSolver(py::module& m)
{
    auto solve1 = static_cast<EquilibriumResult(EquilibriumSolver::*)(ChemicalState&, double, double, VectorConstRef)>(&EquilibriumSolver::solve);
    auto solve2 = static_cast<EquilibriumResult(EquilibriumSolver::*)(ChemicalState&, double, double, const double*)>(&EquilibriumSolver::solve);
    auto solve3 = static_cast<EquilibriumResult(EquilibriumSolver::*)(ChemicalState&, const EquilibriumProblem&)>(&EquilibriumSolver::solve);
    auto solve4 = static_cast<EquilibriumResult(EquilibriumSolver::*)(ChemicalState&)>(&EquilibriumSolver::solve);

    auto approximate1 = static_cast<EquilibriumResult(EquilibriumSolver::*)(ChemicalState&, double, double, VectorConstRef)>(&EquilibriumSolver::approximate);
    auto approximate2 = static_cast<EquilibriumResult(EquilibriumSolver::*)(ChemicalState&, const EquilibriumProblem&)>(&EquilibriumSolver::approximate);
    auto approximate3 = static_cast<EquilibriumResult(EquilibriumSolver::*)(ChemicalState&)>(&EquilibriumSolver::approximate);

    py::class_<EquilibriumSolver>(m, "EquilibriumSolver")
        .def(py::init<const ChemicalSystem&>())
        .def(py::init<const Partition&>())
        .def("setOptions", &EquilibriumSolver::setOptions)
        .def("approximate", approximate1)
        .def("approximate", approximate2)
        .def("approximate", approximate3)
        .def("solve", solve1)
        .def("solve", solve2)
        .def("solve", solve3)
        .def("solve", solve4)
        .def("properties", &EquilibriumSolver::properties, py::return_value_policy::reference_internal)
        .def("sensitivity", &EquilibriumSolver::sensitivity, py::return_value_policy::reference_internal)
        .def("result", &EquilibriumSolver::result, py::return_value_policy::reference_internal)

        // DEPRECATED METHODS: TO BE REMOVED
        .def("setPartition", &EquilibriumSolver::setPartition)
        ;
}

} // namespace Reaktoro
