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

// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

// Reaktoro includes

#include <Reaktoro/Oilphase/oilphase.cpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/GaseousMixture.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/GeneralMixture.hpp>
#include "W:\release\Projects\Reaktoro\demos\cpp\phaseid.hpp"

namespace Reaktoro {

	void exportOil(py::module& m)
	{
		m.def("convertPhase", &convertPhase<OilPhase>);
		m.def("convertPhase", &convertPhase<GaseousPhase>);
		m.def("convertPhase", &convertPhase<AqueousPhase>);

		py::class_<OilSpecies, Species>(m, "OilSpecies")
			.def(py::init<>())
			.def(py::init<const GaseousSpecies&>())
			.def("setCriticalTemperature", &OilSpecies::setCriticalTemperature)
			.def("setCriticalPressure", &OilSpecies::setCriticalPressure)
			.def("setAcentricFactor", &OilSpecies::setAcentricFactor)
			.def("setThermoData", &OilSpecies::setThermoData)
			.def("criticalTemperature", &OilSpecies::criticalTemperature)
			.def("criticalPressure", &OilSpecies::criticalPressure)
			.def("acentricFactor", &OilSpecies::acentricFactor)
			.def("thermoData", &OilSpecies::thermoData, py::return_value_policy::reference_internal)
			;
		
		py::class_<GeneralMixture<OilSpecies>>(m, "_GeneralMixture_oil")
			.def(py::init<>())
			.def(py::init<const std::vector<OilSpecies>&>())
			;

		py::class_<GeneralMixture<GaseousSpecies>>(m, "_GeneralMixture_gas")
			.def(py::init<>())
			.def(py::init<const std::vector<GaseousSpecies>&>())
			;

		py::class_<OilMixture, GeneralMixture<OilSpecies>>(m, "OilMixture")
			.def(py::init<>())
			.def(py::init<const std::vector<OilSpecies>>())
			;
		
		py::class_<GaseousMixture, GeneralMixture<GaseousSpecies>>(m, "GaseousMixture")
			.def(py::init<>())
			.def(py::init<const std::vector<GaseousSpecies>>())
			;

		py::class_<OilPhase, Phase>(m, "OilPhase")
			.def(py::init<>())
			.def(py::init<const OilMixture&>())
			.def("setChemicalModelRedlichKwong", &OilPhase::setChemicalModelRedlichKwong, py::return_value_policy::reference_internal)
			.def("setChemicalModelSoaveRedlichKwong", &OilPhase::setChemicalModelSoaveRedlichKwong, py::return_value_policy::reference_internal)
			.def("setChemicalModelPengRobinson", &OilPhase::setChemicalModelPengRobinson, py::return_value_policy::reference_internal)			
			.def("mixture", &OilPhase::mixture, py::return_value_policy::reference_internal)
			;
			
		py::enum_<phaseIdentificationMethod>(m, "phaseIdentificationMethod")
			.value("VolumeMethod", phaseIdentificationMethod::VolumeMethod)
			.value("CriticalPointMethods", phaseIdentificationMethod::CriticalPointMethods)
			.value("IsothermalCompressibilityMethods", phaseIdentificationMethod::IsothermalCompressibilityMethods)
			.value("workanalisys", phaseIdentificationMethod::workanalisys)
			.value("workanalisysPengRobinson", phaseIdentificationMethod::workanalisysPengRobinson)
			.value("Gibbs_residual_based", phaseIdentificationMethod::Gibbs_residual_based)
			.value("NoMethod", phaseIdentificationMethod::NoMethod)
			;
	}

} // namespace Reaktoro
