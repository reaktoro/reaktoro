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
