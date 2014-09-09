/*
 * Reaktor is a C++ library for computational reaction modelling.
 *
 * Copyright (C) 2014 Allan Leal
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

// Reaktor includes
#include <Reaktor/Thermodynamics/WaterThermoModel.hpp>

namespace Reaktor {

struct WaterHelmholtzState
{
	WaterHelmholtzState();

	/// The specific Helmholtz free energy of water (unit: J/kg)
	double helmholtz;

	/// The first-order partial derivative of the specific Helmholtz free energy of water with respect to temperature
	double helmholtzT;

	/// The first-order partial derivative of the specific Helmholtz free energy of water with respect to density
	double helmholtzD;

	/// The second-order partial derivative of the specific Helmholtz free energy of water with respect to temperature
	double helmholtzTT;

	/// The second-order partial derivative of the specific Helmholtz free energy of water with respect to temperature and density
	double helmholtzTD;

	/// The second-order partial derivative of the specific Helmholtz free energy of water with respect to density
	double helmholtzDD;

	/// The third-order partial derivative of the specific Helmholtz free energy of water with respect to temperature
	double helmholtzTTT;

	/// The third-order partial derivative of the specific Helmholtz free energy of water with respect to temperature, temperature, and density
	double helmholtzTTD;

	/// The third-order partial derivative of the specific Helmholtz free energy of water with respect to temperature, density, and density
	double helmholtzTDD;

	/// The third-order partial derivative of the specific Helmholtz free energy of water with respect to density
	double helmholtzDDD;
};

/**
 * Calculates the Helmholtz free energy state of water with the Wagner and Pruss (1995) equations of state
 * @param T The temperature of water (in units of K)
 * @param D The density of water (in units of kg/m3)
 * @return The Helmholtz free energy state of water
 * @see WaterHelmholtzState
 */
auto waterHelmholtz(double T, double D) -> WaterHelmholtzState;

/**
 * Calculates the Helmholtz free energy state of water
 * @param T The temperature of water (in units of K)
 * @param D The density of water (in units of kg/m3)
 * @param model The thermodynamic model used to calculate the Helmholtz state of water
 * @return The Helmholtz free energy state of water
 * @see WaterHelmholtzState, WaterThermoModel
 */
auto waterHelmholtz(double T, double D, WaterThermoModel model) -> WaterHelmholtzState;

} // namespace Reaktor
