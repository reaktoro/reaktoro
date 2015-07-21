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

#include "PyAqueousPhase.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Mixtures/AqueousMixture.hpp>
#include <Reaktoro/Thermodynamics/Phases/AqueousPhase.hpp>
#include <Reaktoro/Thermodynamics/Species/AqueousSpecies.hpp>

namespace Reaktoro {

auto export_AqueousPhase() -> void
{
    py::class_<AqueousPhase, py::bases<Phase>>("AqueousPhase")
        .def(py::init<>())
        .def(py::init<const AqueousMixture&>())
        .def("setChemicalModelHKF", &AqueousPhase::setChemicalModelHKF)
        .def("setChemicalModelPitzerHMW", &AqueousPhase::setChemicalModelPitzerHMW)
        .def("setActivityModel", &AqueousPhase::setActivityModel)
        .def("setActivityModelIdeal", &AqueousPhase::setActivityModelIdeal)
        .def("setActivityModelSetschenow", &AqueousPhase::setActivityModelSetschenow)
        .def("setActivityModelDuanSunCO2", &AqueousPhase::setActivityModelDuanSunCO2)
        .def("setActivityModelDrummondCO2", &AqueousPhase::setActivityModelDrummondCO2)
        .def("setActivityModelRumpfCO2", &AqueousPhase::setActivityModelRumpfCO2)
        .def("mixture", &AqueousPhase::mixture, py::return_internal_reference<>())
        ;
}

} // namespace Reaktoro
