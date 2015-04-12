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

#pragma once

namespace Reaktoro {

struct WaterThermoState
{
	/// The temperature of water (in units of K)
	double temperature = 0.0;

	/// The specific volume of water (in units of m3/kg)
	double volume = 0.0;

	/// The specific entropy of water (in units of J/(kg*K))
	double entropy = 0.0;

	/// The specific Helmholtz free energy of water (in units of J/kg)
	double helmholtz = 0.0;

	/// The specific internal energy of water (in units of J/kg)
	double internal_energy = 0.0;

	/// The specific enthalpy of water (in units of J/kg)
	double enthalpy = 0.0;

	/// The specific Gibbs free energy of water (in units of J/kg)
	double gibbs = 0.0;

	/// The specific isochoric heat capacity of water (in units of J/(kg*K))
	double cv = 0.0;

	/// The specific isobaric heat capacity of water (in units of J/(kg*K))
	double cp = 0.0;

	/// The specific density of water (in units of kg/m3)
	double density = 0.0;

	/// The first-order partial derivative of density with respect to temperature (in units of (kg/m3)/K)
	double densityT = 0.0;

	/// The first-order partial derivative of density with respect to pressure (in units of (kg/m3)/Pa)
	double densityP = 0.0;

	/// The second-order partial derivative of density with respect to temperature (in units of (kg/m3)/(K*K))
	double densityTT = 0.0;

	/// The second-order partial derivative of density with respect to temperature and pressure (in units of (kg/m3)/(K*Pa))
	double densityTP = 0.0;

	/// The second-order partial derivative of density with respect to pressure (in units of (kg/m3)/(Pa*Pa))
	double densityPP = 0.0;

	/// The pressure of water (in units of Pa)
	double pressure = 0.0;

	/// The first-order partial derivative of pressure with respect to temperature (in units of Pa/K)
	double pressureT = 0.0;

	/// The first-order partial derivative of pressure with respect to density (in units of Pa/(kg/m3))
	double pressureD = 0.0;

	/// The second-order partial derivative of pressure with respect to temperature (in units of Pa/(K*K))
	double pressureTT = 0.0;

	/// The second-order partial derivative of pressure with respect to temperature and density (in units of Pa/(K*kg/m3))
	double pressureTD = 0.0;

	/// The second-order partial derivative of pressure with respect to density (in units of Pa/((kg/m3)*(kg/m3)))
	double pressureDD = 0.0;
};

} // namespace Reaktoro
