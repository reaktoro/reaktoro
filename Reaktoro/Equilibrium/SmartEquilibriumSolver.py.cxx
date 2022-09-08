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
#include <Reaktoro/Equilibrium/EquilibriumRestrictions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumSolver.hpp>
using namespace Reaktoro;

void exportSmartEquilibriumSolver(py::module& m)
{
    py::class_<SmartEquilibriumSolver>(m, "SmartEquilibriumSolver")
        .def(py::init<ChemicalSystem const&>())
        .def(py::init<EquilibriumSpecs const&>())

        .def("solve", py::overload_cast<ChemicalState&>(&SmartEquilibriumSolver::solve), "Equilibrate a chemical state.", py::arg("state"))
        .def("solve", py::overload_cast<ChemicalState&, EquilibriumRestrictions const&>(&SmartEquilibriumSolver::solve), "Equilibrate a chemical state respecting given reactivity restrictions.", py::arg("state"), py::arg("restrictions"))
        .def("solve", py::overload_cast<ChemicalState&, EquilibriumConditions const&>(&SmartEquilibriumSolver::solve), "Equilibrate a chemical state respecting given constraint conditions.", py::arg("state"), py::arg("conditions"))
        .def("solve", py::overload_cast<ChemicalState&, EquilibriumConditions const&, EquilibriumRestrictions const&>(&SmartEquilibriumSolver::solve), "Equilibrate a chemical state respecting given constraint conditions and reactivity restrictions.", py::arg("state"), py::arg("conditions"), py::arg("restrictions"))

        .def("solve", py::overload_cast<ChemicalState&, EquilibriumSensitivity&>(&SmartEquilibriumSolver::solve), "Equilibrate a chemical state and compute sensitivity derivatives.", py::arg("state"), py::arg("sensitivity"))
        .def("solve", py::overload_cast<ChemicalState&, EquilibriumSensitivity&, EquilibriumRestrictions const&>(&SmartEquilibriumSolver::solve), "Equilibrate a chemical state respecting given reactivity restrictions and compute sensitivity derivatives.", py::arg("state"), py::arg("sensitivity"), py::arg("restrictions"))
        .def("solve", py::overload_cast<ChemicalState&, EquilibriumSensitivity&, EquilibriumConditions const&>(&SmartEquilibriumSolver::solve), "Equilibrate a chemical state respecting given constraint conditions and compute sensitivity derivatives.", py::arg("state"), py::arg("sensitivity"), py::arg("conditions"))
        .def("solve", py::overload_cast<ChemicalState&, EquilibriumSensitivity&, EquilibriumConditions const&, EquilibriumRestrictions const&>(&SmartEquilibriumSolver::solve), "Equilibrate a chemical state respecting given constraint conditions and reactivity restrictions and compute sensitivity derivatives.", py::arg("state"), py::arg("sensitivity"), py::arg("conditions"), py::arg("restrictions"))

        .def("setOptions", &SmartEquilibriumSolver::setOptions)
        ;
}
