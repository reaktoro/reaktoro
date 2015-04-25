// Reaktoro is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
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
    py::class_<GaseousPhase>("GaseousPhase")
        .def(py::init<>())
        .def(py::init<const std::vector<GaseousSpecies>&>())
        .def("setActivityModel", &GaseousPhase::setActivityModel)
        .def("setActivityModelIdeal", &GaseousPhase::setActivityModelIdeal)
        .def("setActivityModelDuanSunCO2", &GaseousPhase::setActivityModelDuanSunCO2)
        .def("setActivityModelSpycherPruessH2OCO2", &GaseousPhase::setActivityModelSpycherPruessH2OCO2)
        .def("setActivityModelSpycherReedH2OCO2CH4", &GaseousPhase::setActivityModelSpycherReedH2OCO2CH4)
        .def("setActivityModelPengRobinson", &GaseousPhase::setActivityModelPengRobinson)
        .def("concentrations", &GaseousPhase::concentrations)
        .def("activityCoefficients", &GaseousPhase::activityCoefficients)
        .def("activities", &GaseousPhase::activities)
        ;
}

} // namespace Reaktoro
