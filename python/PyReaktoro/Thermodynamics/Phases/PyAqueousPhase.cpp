// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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
#include <Reaktoro/Extensions/Geochemistry/AqueousMixture.hpp>
#include <Reaktoro/Models/ActivityModelDebyeHuckel.hpp>
#include <Reaktoro/Thermodynamics/Phases/AqueousPhase.hpp>
#include <Reaktoro/Thermodynamics/Species/AqueousSpecies.hpp>

namespace Reaktoro {

void exportAqueousPhase(py::module& m)
{
	auto setChemicalModelDebyeHuckel1 = static_cast<AqueousPhase&(AqueousPhase::*)()>(&AqueousPhase::setChemicalModelDebyeHuckel);
	auto setChemicalModelDebyeHuckel2 = static_cast<AqueousPhase&(AqueousPhase::*)(const DebyeHuckelParams&)>(&AqueousPhase::setChemicalModelDebyeHuckel);

    py::class_<AqueousPhase, Phase>(m, "AqueousPhase")
        .def(py::init<>())
        .def(py::init<const AqueousMixture&>())
        .def("setInterpolationPoints", &AqueousPhase::setInterpolationPoints, py::return_value_policy::reference_internal)
        .def("setChemicalModelIdeal", &AqueousPhase::setChemicalModelIdeal, py::return_value_policy::reference_internal)
        .def("setChemicalModelDebyeHuckel", setChemicalModelDebyeHuckel1, py::return_value_policy::reference_internal)
        .def("setChemicalModelDebyeHuckel", setChemicalModelDebyeHuckel2, py::return_value_policy::reference_internal)
        .def("setChemicalModelHKF", &AqueousPhase::setChemicalModelHKF, py::return_value_policy::reference_internal)
        .def("setChemicalModelPitzerHMW", &AqueousPhase::setChemicalModelPitzerHMW, py::return_value_policy::reference_internal)
        .def("setActivityModel", &AqueousPhase::setActivityModel, py::return_value_policy::reference_internal)
        .def("setActivityModelIdeal", &AqueousPhase::setActivityModelIdeal, py::return_value_policy::reference_internal)
        .def("setActivityModelSetschenow", &AqueousPhase::setActivityModelSetschenow, py::return_value_policy::reference_internal)
        .def("setActivityModelDuanSunCO2", &AqueousPhase::setActivityModelDuanSunCO2, py::return_value_policy::reference_internal)
        .def("setActivityModelDrummondCO2", &AqueousPhase::setActivityModelDrummondCO2, py::return_value_policy::reference_internal)
        .def("setActivityModelRumpfCO2", &AqueousPhase::setActivityModelRumpfCO2, py::return_value_policy::reference_internal)
        .def("mixture", &AqueousPhase::mixture, py::return_value_policy::reference_internal)
        ;
}

} // namespace Reaktoro
