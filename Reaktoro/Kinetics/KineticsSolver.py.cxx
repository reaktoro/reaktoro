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
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConditions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumRestrictions.hpp>
#include <Reaktoro/Kinetics/KineticsOptions.hpp>
#include <Reaktoro/Kinetics/KineticsResult.hpp>
#include <Reaktoro/Kinetics/KineticsSensitivity.hpp>
#include <Reaktoro/Kinetics/KineticsSolver.hpp>
using namespace Reaktoro;

void exportKineticsSolver(py::module& m)
{
    py::class_<KineticsSolver>(m, "KineticsSolver")
        .def(py::init<ChemicalSystem const&>())
        .def(py::init<EquilibriumSpecs const&>())

        // .def("solve", py::overload_cast<ChemicalState&>(&KineticsSolver::solve))
        // .def("solve", py::overload_cast<ChemicalState&, EquilibriumRestrictions const&>(&KineticsSolver::solve))
        // .def("solve", py::overload_cast<ChemicalState&, EquilibriumConditions const&>(&KineticsSolver::solve))
        // .def("solve", py::overload_cast<ChemicalState&, EquilibriumConditions const&, EquilibriumRestrictions const&>(&KineticsSolver::solve))

        // .def("solve", py::overload_cast<ChemicalState&, EquilibriumSensitivity&>(&KineticsSolver::solve))
        // .def("solve", py::overload_cast<ChemicalState&, EquilibriumSensitivity&, EquilibriumRestrictions const&>(&KineticsSolver::solve))
        // .def("solve", py::overload_cast<ChemicalState&, EquilibriumSensitivity&, EquilibriumConditions const&>(&KineticsSolver::solve))
        // .def("solve", py::overload_cast<ChemicalState&, EquilibriumSensitivity&, EquilibriumConditions const&, EquilibriumRestrictions const&>(&KineticsSolver::solve))

        // .def("solve", py::overload_cast<ChemicalState&, ArrayXdConstRef>(&KineticsSolver::solve))
        // .def("solve", py::overload_cast<ChemicalState&, EquilibriumRestrictions const&, ArrayXdConstRef>(&KineticsSolver::solve))
        // .def("solve", py::overload_cast<ChemicalState&, EquilibriumConditions const&, ArrayXdConstRef>(&KineticsSolver::solve))
        // .def("solve", py::overload_cast<ChemicalState&, EquilibriumConditions const&, EquilibriumRestrictions const&, ArrayXdConstRef>(&KineticsSolver::solve))

        .def("setOptions", &KineticsSolver::setOptions)
        ;
}
