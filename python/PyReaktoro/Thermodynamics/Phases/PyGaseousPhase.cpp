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

#include "PyGaseousPhase.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Mixtures/GaseousMixture.hpp>
#include <Reaktoro/Thermodynamics/Phases/GaseousPhase.hpp>
#include <Reaktoro/Thermodynamics/Species/GaseousSpecies.hpp>

namespace Reaktoro {

auto export_GaseousPhase() -> void
{
    py::class_<GaseousPhase, py::bases<Phase>>("GaseousPhase")
        .def(py::init<>())
        .def(py::init<const GaseousMixture&>())
        .def("setChemicalModelIdeal", &GaseousPhase::setChemicalModelIdeal, py::return_internal_reference<>())
        .def("setChemicalModelVanDerWaals", &GaseousPhase::setChemicalModelVanDerWaals, py::return_internal_reference<>())
        .def("setChemicalModelRedlichKwong", &GaseousPhase::setChemicalModelRedlichKwong, py::return_internal_reference<>())
        .def("setChemicalModelSoaveRedlichKwong", &GaseousPhase::setChemicalModelSoaveRedlichKwong, py::return_internal_reference<>())
        .def("setChemicalModelPengRobinson", &GaseousPhase::setChemicalModelPengRobinson, py::return_internal_reference<>())
        .def("setChemicalModelSpycherPruessEnnis", &GaseousPhase::setChemicalModelSpycherPruessEnnis, py::return_internal_reference<>())
        .def("setChemicalModelSpycherReed", &GaseousPhase::setChemicalModelSpycherReed, py::return_internal_reference<>())
        .def("mixture", &GaseousPhase::mixture, py::return_internal_reference<>())
        ;
}

} // namespace Reaktoro
