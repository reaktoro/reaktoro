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

#include "HydrocarbonSpecies.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>

namespace Reaktoro {

struct HydrocarbonSpecies::Impl
{
	// The critical temperature of the hydrocarbon species (in units of K)
	double critical_temperature = 0.0;

	// The critical pressure of the hydrocarbon species (in units of Pa)
	double critical_pressure = 0.0;

	// The acentric factor of the hydrocarbon species
	double acentric_factor = 0.0;

	/// The termidynamic data of the hydrocarbon species
	HydrocarbonSpeciesThermoData thermo;
};

HydrocarbonSpecies::HydrocarbonSpecies()
	: pimpl(new Impl())
{}

HydrocarbonSpecies::HydrocarbonSpecies(const Species& species)
	: Species(species), pimpl(new Impl())
{}

auto HydrocarbonSpecies::setCriticalTemperature(double val) -> void
{
	Assert(val > 0.0, "Cannot set the critical temperature of the hydrocarbon specie `" + name() + "`.",
		"The given critical temperature `" + std::to_string(val) + "` is not positive.");
	pimpl->critical_temperature = val;
}

auto HydrocarbonSpecies::setCriticalPressure(double val) -> void
{
	Assert(val > 0.0, "Cannot set the critical pressure of the hydrocarbon specie `" + name() + "`.",
		"The given critical pressure `" + std::to_string(val) + "` is not positive.");
	pimpl->critical_pressure = val;
}

auto HydrocarbonSpecies::setAcentricFactor(double val) -> void
{
	pimpl->acentric_factor = val;
}

auto HydrocarbonSpecies::setThermoData(const HydrocarbonSpeciesThermoData& thermo) -> void
{
	pimpl->thermo = thermo;
}

auto HydrocarbonSpecies::criticalTemperature() const -> double
{
	return pimpl->critical_temperature;
}

auto HydrocarbonSpecies::criticalPressure() const -> double
{
	return pimpl->critical_pressure;
}

auto HydrocarbonSpecies::acentricFactor() const -> double
{
	return pimpl->acentric_factor;
}

auto HydrocarbonSpecies::thermoData() const -> const HydrocarbonSpeciesThermoData&
{
	return pimpl->thermo;
}



} // namespace Reaktoro
