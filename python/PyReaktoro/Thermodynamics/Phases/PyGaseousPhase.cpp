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

// pybind11 includes
#include <pybind11/pybind11.h>
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Mixtures/GaseousMixture.hpp>
#include <Reaktoro/Thermodynamics/Phases/GaseousPhase.hpp>
#include <Reaktoro/Thermodynamics/Species/GaseousSpecies.hpp>

namespace Reaktoro {

void exportGaseousPhase(py::module& m)
{
    py::class_<GaseousPhase, Phase>(m, "GaseousPhase")
        .def(py::init<>())
        .def(py::init<const GaseousMixture&>())
        .def("setChemicalModelIdeal", &GaseousPhase::setChemicalModelIdeal, py::return_value_policy::reference_internal)
        .def("setChemicalModelVanDerWaals", &GaseousPhase::setChemicalModelVanDerWaals, py::return_value_policy::reference_internal)
        .def("setChemicalModelRedlichKwong", &GaseousPhase::setChemicalModelRedlichKwong, py::return_value_policy::reference_internal)
        .def("setChemicalModelSoaveRedlichKwong", &GaseousPhase::setChemicalModelSoaveRedlichKwong, py::return_value_policy::reference_internal)
        .def("setChemicalModelPengRobinson", &GaseousPhase::setChemicalModelPengRobinson, py::return_value_policy::reference_internal)
        .def("setChemicalModelSpycherPruessEnnis", &GaseousPhase::setChemicalModelSpycherPruessEnnis, py::return_value_policy::reference_internal)
        .def("setChemicalModelSpycherReed", &GaseousPhase::setChemicalModelSpycherReed, py::return_value_policy::reference_internal)
        .def("mixture", &GaseousPhase::mixture, py::return_value_policy::reference_internal)
        ;
}

} // namespace Reaktoro
