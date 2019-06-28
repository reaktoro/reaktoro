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

#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/GaseousMixture.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/GeneralMixture.hpp>
#include "W:\release\Projects\Reaktoro\demos\cpp\phaseid.hpp"
#include <reaktoro/Thermodynamics/Mixtures/LiquidMixture.hpp>

namespace Reaktoro {

	void exportOil(py::module& m)
	{
		py::class_<GeneralMixture<LiquidSpecies>>(m, "_GeneralMixture_oil")
			.def(py::init<>())
			.def(py::init<const std::vector<LiquidSpecies>&>())
			;

		py::class_<GeneralMixture<GaseousSpecies>>(m, "_GeneralMixture_gas")
			.def(py::init<>())
			.def(py::init<const std::vector<GaseousSpecies>&>())
			;

		py::class_<LiquidMixture, GeneralMixture<LiquidSpecies>>(m, "OilMixture")
			.def(py::init<>())
			.def(py::init<const std::vector<LiquidSpecies>>())
			;
		
		py::class_<GaseousMixture, GeneralMixture<GaseousSpecies>>(m, "GaseousMixture")
			.def(py::init<>())
			.def(py::init<const std::vector<GaseousSpecies>>())
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
