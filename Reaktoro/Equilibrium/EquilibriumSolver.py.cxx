// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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
#include <Reaktoro/pybind11.hxx>

// Reaktoro includes
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConditions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumRestrictions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>
using namespace Reaktoro;

void exportEquilibriumSolver(py::module& m)
{
    py::class_<EquilibriumSolver>(m, "EquilibriumSolver")
        .def(py::init<const ChemicalSystem&>())
        .def(py::init<const EquilibriumSpecs&>())

        .def("solve", py::overload_cast<ChemicalState&>(&EquilibriumSolver::solve))
        .def("solve", py::overload_cast<ChemicalState&, const EquilibriumRestrictions&>(&EquilibriumSolver::solve))
        .def("solve", py::overload_cast<ChemicalState&, const EquilibriumConditions&>(&EquilibriumSolver::solve))
        .def("solve", py::overload_cast<ChemicalState&, const EquilibriumConditions&, const EquilibriumRestrictions&>(&EquilibriumSolver::solve))

        // .def("solve", py::overload_cast<ChemicalState&, EquilibriumSensitivity&>(&EquilibriumSolver::solve))
        // .def("solve", py::overload_cast<ChemicalState&, EquilibriumSensitivity&, const EquilibriumRestrictions&>(&EquilibriumSolver::solve))
        // .def("solve", py::overload_cast<ChemicalState&, EquilibriumSensitivity&, const EquilibriumConditions&>(&EquilibriumSolver::solve))
        // .def("solve", py::overload_cast<ChemicalState&, EquilibriumSensitivity&, const EquilibriumConditions&, const EquilibriumRestrictions&>(&EquilibriumSolver::solve))

        .def("solve", py::overload_cast<ChemicalState&, ArrayXdConstRef>(&EquilibriumSolver::solve))
        .def("solve", py::overload_cast<ChemicalState&, const EquilibriumRestrictions&, ArrayXdConstRef>(&EquilibriumSolver::solve))
        .def("solve", py::overload_cast<ChemicalState&, const EquilibriumConditions&, ArrayXdConstRef>(&EquilibriumSolver::solve))
        .def("solve", py::overload_cast<ChemicalState&, const EquilibriumConditions&, const EquilibriumRestrictions&, ArrayXdConstRef>(&EquilibriumSolver::solve))

        .def("setOptions", &EquilibriumSolver::setOptions)
        .def("conservativeMatrix", &EquilibriumSolver::conservativeMatrix)
        .def("componentAmounts", &EquilibriumSolver::componentAmounts)
        ;
}
