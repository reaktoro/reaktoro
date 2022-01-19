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
#include <Reaktoro/Equilibrium/EquilibriumProblem.hpp>
using namespace Reaktoro;

void exportEquilibriumProblem(py::module& m)
{
    py::class_<EquilibriumProblem, EquilibriumConditions, EquilibriumRestrictions>(m, "EquilibriumProblem")
        .def(py::init<const ChemicalSystem&>())
        .def(py::init<const EquilibriumSpecs&>())
        .def("startWithTemperature", py::overload_cast<real>(&EquilibriumProblem::startWithTemperature))
        .def("startWithTemperature", py::overload_cast<real, String>(&EquilibriumProblem::startWithTemperature))
        .def("startWithPressure", py::overload_cast<real>(&EquilibriumProblem::startWithPressure))
        .def("startWithPressure", py::overload_cast<real, String>(&EquilibriumProblem::startWithPressure))
        .def("startWithSpeciesAmounts", &EquilibriumProblem::startWithSpeciesAmounts)
        .def("startWith", py::overload_cast<String, real, String>(&EquilibriumProblem::startWith))
        .def("startWith", py::overload_cast<Index, real, String>(&EquilibriumProblem::startWith))
        .def("startWithState", &EquilibriumProblem::startWithState)
        .def("startWithComponentAmounts", &EquilibriumProblem::startWithComponentAmounts)
        .def("initialTemperature", &EquilibriumProblem::initialTemperature)
        .def("initialPressure", &EquilibriumProblem::initialPressure)
        .def("initialSpeciesAmounts", &EquilibriumProblem::initialSpeciesAmounts, return_internal_ref)
        .def("initialComponentAmounts", &EquilibriumProblem::initialComponentAmounts, return_internal_ref)
        ;
}
