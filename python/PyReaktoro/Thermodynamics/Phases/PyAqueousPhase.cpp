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

#include "PyAqueousPhase.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Mixtures/AqueousMixture.hpp>
#include <Reaktoro/Thermodynamics/Models/AqueousChemicalModelDebyeHuckel.hpp>
#include <Reaktoro/Thermodynamics/Phases/AqueousPhase.hpp>
#include <Reaktoro/Thermodynamics/Species/AqueousSpecies.hpp>

namespace Reaktoro {

auto export_AqueousPhase() -> void
{
	auto setChemicalModelDebyeHuckel1 = static_cast<AqueousPhase&(AqueousPhase::*)()>(&AqueousPhase::setChemicalModelDebyeHuckel);
	auto setChemicalModelDebyeHuckel2 = static_cast<AqueousPhase&(AqueousPhase::*)(const DebyeHuckelParams&)>(&AqueousPhase::setChemicalModelDebyeHuckel);

    py::class_<AqueousPhase, py::bases<Phase>>("AqueousPhase")
        .def(py::init<>())
        .def(py::init<const AqueousMixture&>())
        .def("setInterpolationPoints", &AqueousPhase::setInterpolationPoints, py::return_internal_reference<>())
        .def("setChemicalModelIdeal", &AqueousPhase::setChemicalModelIdeal, py::return_internal_reference<>())
        .def("setChemicalModelDebyeHuckel", setChemicalModelDebyeHuckel1, py::return_internal_reference<>())
        .def("setChemicalModelDebyeHuckel", setChemicalModelDebyeHuckel2, py::return_internal_reference<>())
        .def("setChemicalModelHKF", &AqueousPhase::setChemicalModelHKF, py::return_internal_reference<>())
        .def("setChemicalModelPitzerHMW", &AqueousPhase::setChemicalModelPitzerHMW, py::return_internal_reference<>())
        .def("setActivityModel", &AqueousPhase::setActivityModel, py::return_internal_reference<>())
        .def("setActivityModelIdeal", &AqueousPhase::setActivityModelIdeal, py::return_internal_reference<>())
        .def("setActivityModelSetschenow", &AqueousPhase::setActivityModelSetschenow, py::return_internal_reference<>())
        .def("setActivityModelDuanSunCO2", &AqueousPhase::setActivityModelDuanSunCO2, py::return_internal_reference<>())
        .def("setActivityModelDrummondCO2", &AqueousPhase::setActivityModelDrummondCO2, py::return_internal_reference<>())
        .def("setActivityModelRumpfCO2", &AqueousPhase::setActivityModelRumpfCO2, py::return_internal_reference<>())
        .def("mixture", &AqueousPhase::mixture, py::return_internal_reference<>())
        ;
}

} // namespace Reaktoro
