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
#include <Reaktoro/Thermodynamics/Species/LiquidSpecies.hpp>

namespace Reaktoro {

	void exportLiquidSpecies(py::module& m)
	{
		py::class_<LiquidSpecies, Species>(m, "LiquidSpecies")
			.def(py::init<>())
			.def(py::init<const GaseousSpecies&>())
			.def("setCriticalTemperature", &LiquidSpecies::setCriticalTemperature)
			.def("setCriticalPressure", &LiquidSpecies::setCriticalPressure)
			.def("setAcentricFactor", &LiquidSpecies::setAcentricFactor)
			.def("setThermoData", &LiquidSpecies::setThermoData)
			.def("criticalTemperature", &LiquidSpecies::criticalTemperature)
			.def("criticalPressure", &LiquidSpecies::criticalPressure)
			.def("acentricFactor", &LiquidSpecies::acentricFactor)
			.def("thermoData", &LiquidSpecies::thermoData, py::return_value_policy::reference_internal)
			;
	}

} // namespace Reaktoro
