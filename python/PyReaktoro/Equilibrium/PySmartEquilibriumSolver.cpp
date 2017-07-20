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

#include "PyEquilibriumSolver.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumSolver.hpp>

namespace Reaktoro {

auto export_SmartEquilibriumSolver() -> void
{
    auto learn1 = static_cast<EquilibriumResult(SmartEquilibriumSolver::*)(ChemicalState&, double, double, const Vector&)>(&SmartEquilibriumSolver::learn);
    auto learn2 = static_cast<EquilibriumResult(SmartEquilibriumSolver::*)(ChemicalState&, const EquilibriumProblem&)>(&SmartEquilibriumSolver::learn);

    auto estimate1 = static_cast<EquilibriumResult(SmartEquilibriumSolver::*)(ChemicalState&, double, double, const Vector&)>(&SmartEquilibriumSolver::estimate);
    auto estimate2 = static_cast<EquilibriumResult(SmartEquilibriumSolver::*)(ChemicalState&, const EquilibriumProblem&)>(&SmartEquilibriumSolver::estimate);

    auto solve1 = static_cast<EquilibriumResult(SmartEquilibriumSolver::*)(ChemicalState&, double, double, const Vector&)>(&SmartEquilibriumSolver::solve);
    auto solve2 = static_cast<EquilibriumResult(SmartEquilibriumSolver::*)(ChemicalState&, const EquilibriumProblem&)>(&SmartEquilibriumSolver::solve);

    py::class_<SmartEquilibriumSolver>("SmartEquilibriumSolver", py::no_init)
        .def(py::init<const ChemicalSystem&>())
        .def("setOptions", &SmartEquilibriumSolver::setOptions)
        .def("setPartition", &SmartEquilibriumSolver::setPartition)
        .def("learn", learn1)
        .def("learn", learn2)
        .def("estimate", estimate1)
        .def("estimate", estimate2)
        .def("solve", solve1)
        .def("solve", solve2)
        ;
}

} // namespace Reaktoro
