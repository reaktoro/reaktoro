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
#include <Reaktoro/Equilibrium/SmartEquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumSolver.hpp>

namespace Reaktoro {

void exportSmartEquilibriumSolver(py::module& m)
{
    auto solve1 = static_cast<SmartEquilibriumResult(SmartEquilibriumSolver::*)(ChemicalState&, double, double, VectorConstRef)>(&SmartEquilibriumSolver::solve);
    auto solve2 = static_cast<SmartEquilibriumResult(SmartEquilibriumSolver::*)(ChemicalState&, const EquilibriumProblem&)>(&SmartEquilibriumSolver::solve);

    py::class_<SmartEquilibriumSolver>(m, "SmartEquilibriumSolver")
        .def(py::init<const ChemicalSystem&>())
        .def(py::init<const Partition&>())
        .def("setOptions", &SmartEquilibriumSolver::setOptions)
        .def("solve", solve1)
        .def("solve", solve2)
        .def("properties", &SmartEquilibriumSolver::properties, py::return_value_policy::reference_internal)
        .def("result", &SmartEquilibriumSolver::result, py::return_value_policy::reference_internal)
        .def("outputClusterInfo", &SmartEquilibriumSolver::outputClusterInfo)

        ;
}

} // namespace Reaktoro
