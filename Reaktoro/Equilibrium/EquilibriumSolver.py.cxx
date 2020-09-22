// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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

// pybind11 includes
#include <pybind11/pybind11.h>
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConstraints.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
using namespace Reaktoro;

void exportEquilibriumSolver(py::module& m)
{
    py::class_<EquilibriumSolver>(m, "EquilibriumSolver")
        .def(py::init<const ChemicalSystem&>())
        .def(py::init<const EquilibriumConstraints&>())
        .def("setOptions", &EquilibriumSolver::setOptions)
        .def("constraints", py::overload_cast<>(&EquilibriumSolver::constraints, py::const_))
        .def("constraints", py::overload_cast<>(&EquilibriumSolver::constraints))
        .def("solve", &EquilibriumSolver::solve)
        ;
}
