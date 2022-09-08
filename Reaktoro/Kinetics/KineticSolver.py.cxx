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
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConditions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumRestrictions.hpp>
#include <Reaktoro/Kinetics/KineticSensitivity.hpp>
#include <Reaktoro/Kinetics/KineticSolver.hpp>
using namespace Reaktoro;

void exportKineticSolver(py::module& m)
{
    py::class_<KineticSolver>(m, "KineticSolver")
        .def(py::init<ChemicalSystem const&>())
        .def(py::init<EquilibriumSpecs const&>())

        // .def("solve", py::overload_cast<ChemicalState&>(&KineticSolver::solve))
        // .def("solve", py::overload_cast<ChemicalState&, EquilibriumRestrictions const&>(&KineticSolver::solve))
        // .def("solve", py::overload_cast<ChemicalState&, EquilibriumConditions const&>(&KineticSolver::solve))
        // .def("solve", py::overload_cast<ChemicalState&, EquilibriumConditions const&, EquilibriumRestrictions const&>(&KineticSolver::solve))

        // .def("solve", py::overload_cast<ChemicalState&, EquilibriumSensitivity&>(&KineticSolver::solve))
        // .def("solve", py::overload_cast<ChemicalState&, EquilibriumSensitivity&, EquilibriumRestrictions const&>(&KineticSolver::solve))
        // .def("solve", py::overload_cast<ChemicalState&, EquilibriumSensitivity&, EquilibriumConditions const&>(&KineticSolver::solve))
        // .def("solve", py::overload_cast<ChemicalState&, EquilibriumSensitivity&, EquilibriumConditions const&, EquilibriumRestrictions const&>(&KineticSolver::solve))

        // .def("solve", py::overload_cast<ChemicalState&, ArrayXdConstRef>(&KineticSolver::solve))
        // .def("solve", py::overload_cast<ChemicalState&, EquilibriumRestrictions const&, ArrayXdConstRef>(&KineticSolver::solve))
        // .def("solve", py::overload_cast<ChemicalState&, EquilibriumConditions const&, ArrayXdConstRef>(&KineticSolver::solve))
        // .def("solve", py::overload_cast<ChemicalState&, EquilibriumConditions const&, EquilibriumRestrictions const&, ArrayXdConstRef>(&KineticSolver::solve))

        .def("setOptions", &KineticSolver::setOptions)
        ;
}
