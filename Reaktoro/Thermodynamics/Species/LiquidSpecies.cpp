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

#include "LiquidSpecies.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>

namespace Reaktoro {

struct LiquidSpecies::Impl
{
	// The critical temperature of the liquid species (in units of K)
	double critical_temperature = 0.0;

	// The critical pressure of the liquid species (in units of Pa)
	double critical_pressure = 0.0;

	// The acentric factor of the liquid species
	double acentric_factor = 0.0;

	/// The termidynamic data of the liquid species
    FluidSpeciesThermoData thermo;
};

LiquidSpecies::LiquidSpecies()
	: pimpl(new Impl())
{}

LiquidSpecies::LiquidSpecies(const Species& species)
	: Species(species), pimpl(new Impl())
{}

LiquidSpecies::LiquidSpecies(const GaseousSpecies& species)
	: Species(species) , pimpl(new Impl())
{
	this->setCriticalTemperature(species.criticalTemperature());
	this->setCriticalPressure(species.criticalPressure());
	this->setAcentricFactor(species.acentricFactor());

	auto const& gas_thermo_data = species.thermoData();

	Optional<FluidSpeciesThermoParamsHKF> oil_hkf;
	if (!gas_thermo_data.hkf.empty()) {
        FluidSpeciesThermoParamsHKF const& hkf = gas_thermo_data.hkf.get();
		oil_hkf = FluidSpeciesThermoParamsHKF{ hkf.Gf, hkf.Hf, hkf.Sr, hkf.a, hkf.b, hkf.c, hkf.Tmax };
	}

	this->setThermoData(FluidSpeciesThermoData{
		gas_thermo_data.properties,
		gas_thermo_data.reaction,
		oil_hkf,
		gas_thermo_data.phreeqc
		});
}


auto LiquidSpecies::setCriticalTemperature(double val) -> void
{
	Assert(val > 0.0, "Cannot set the critical temperature of the liquid specie `" + name() + "`.",
		"The given critical temperature `" + std::to_string(val) + "` is not positive.");
	pimpl->critical_temperature = val;
}

auto LiquidSpecies::setCriticalPressure(double val) -> void
{
	Assert(val > 0.0, "Cannot set the critical pressure of the liquid specie `" + name() + "`.",
		"The given critical pressure `" + std::to_string(val) + "` is not positive.");
	pimpl->critical_pressure = val;
}

auto LiquidSpecies::setAcentricFactor(double val) -> void
{
	pimpl->acentric_factor = val;
}

auto LiquidSpecies::setThermoData(const FluidSpeciesThermoData& thermo) -> void
{
	pimpl->thermo = thermo;
}

auto LiquidSpecies::criticalTemperature() const -> double
{
	return pimpl->critical_temperature;
}

auto LiquidSpecies::criticalPressure() const -> double
{
	return pimpl->critical_pressure;
}

auto LiquidSpecies::acentricFactor() const -> double
{
	return pimpl->acentric_factor;
}

auto LiquidSpecies::thermoData() const -> const FluidSpeciesThermoData&
{
	return pimpl->thermo;
}



} // namespace Reaktoro
