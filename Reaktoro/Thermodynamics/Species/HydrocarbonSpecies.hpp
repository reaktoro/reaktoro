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

#pragma once

// C++ includes
#include <string>

// Reaktoro includes
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Thermodynamics/Species/ThermoData.hpp>

namespace Reaktoro {

	/// A type to describe the attributes of a hydrocarbon species
	class HydrocarbonSpecies : public Species
	{
	public:
		/// Construct a default HydrocarbonSpecies instance
		HydrocarbonSpecies();

		/// Construct an HydrocarbonSpecies instance from a Species instance
		HydrocarbonSpecies(const Species& species);

		/// Set the critical temperature of the hydrocarbon species (in units of K)
		auto setCriticalTemperature(double val) -> void;

		/// Set the critical pressure of the hydrocarbon species (in units of Pa)
		auto setCriticalPressure(double val) -> void;

		/// Set the acentric factor of the hydrocarbon species
		auto setAcentricFactor(double val) -> void;

		/// Set the thermodynamic data of the hydrocarbon species.
		auto setThermoData(const HydrocarbonSpeciesThermoData& thermo) -> void;

		/// Return the critical temperature of the hydrocarbon species (in units of K)
		auto criticalTemperature() const -> double;

		/// Return the critical pressure of the hydrocarbon species (in units of Pa)
		auto criticalPressure() const -> double;

		/// Return the acentric factor of the hydrocarbon species
		auto acentricFactor() const -> double;

		/// Return the thermodynamic data of the hydrocarbon species.
		auto thermoData() const -> const HydrocarbonSpeciesThermoData&;

	private:
		struct Impl;

		std::shared_ptr<Impl> pimpl;
	};

} // namespace Reaktoro
