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
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>

namespace Reaktoro {

auto export_EquilibriumSolver() -> void
{
    auto setPartition1 = static_cast<void(EquilibriumSolver::*)(const Partition&)>(&EquilibriumSolver::setPartition);
    auto setPartition2 = static_cast<void(EquilibriumSolver::*)(std::string)>(&EquilibriumSolver::setPartition);

    py::class_<EquilibriumSolver>("EquilibriumSolver", py::no_init)
        .def(py::init<const ChemicalSystem&>())
        .def("setOptions", &EquilibriumSolver::setOptions)
        .def("setPartition", setPartition1)
        .def("setPartition", setPartition2)
        .def("approximate", &EquilibriumSolver::approximate)
        .def("solve", &EquilibriumSolver::solve)
        ;
}

} // namespace Reaktoro
