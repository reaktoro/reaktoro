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

// C++ includes
#include <iostream>

// Reaktor includes
#include <Reaktor/Thermo/WaterThermoModel.hpp>

namespace Reaktor {

// Forward declarations
struct WaterHelmholtz;

struct WaterThermo
{
	/// Constructs a default @ref WaterThermo instance
	WaterThermo();

	/// The temperature of water in units of K.
	double temperature;

	/// The specific volume of water in units of m3/kg.
	double volume;

	/// The specific entropy of water in units of J/(kg:K).
	double entropy;

	/// The specific Helmholtz free energy of water in units of J/kg.
	double helmholtz;

	/// The specific internal energy of water in units of J/kg.
	double internal_energy;

	/// The specific enthalpy of water in units of J/kg.
	double enthalpy;

	/// The specific Gibbs free energy of water in units of J/kg.
	double gibbs;

	/// The specific isochoric heat capacity of water in units of J/(kg:K).
	double cv;

	/// The specific isobaric heat capacity of water in units of J/(kg:K).
	double cp;

	/// The specific density of water in units of kg/m3.
	double density;

	/// The first-order partial derivative of density with respect to temperature in units of (kg/m3)/K.
	double densityT;

	/// The first-order partial derivative of density with respect to pressure in units of (kg/m3)/Pa.
	double densityP;

	/// The second-order partial derivative of density with respect to temperature in units of (kg/m3)/(K:K).
	double densityTT;

	/// The second-order partial derivative of density with respect to temperature and pressure in units of (kg/m3)/(K:Pa).
	double densityTP;

	/// The second-order partial derivative of density with respect to pressure in units of (kg/m3)/(Pa:Pa).
	double densityPP;

	/// The pressure of water in units of Pa.
	double pressure;

	/// The first-order partial derivative of pressure with respect to temperature in units of Pa/K.
	double pressureT;

	/// The first-order partial derivative of pressure with respect to density in units of Pa/(kg/m3).
	double pressureD;

	/// The second-order partial derivative of pressure with respect to temperature in units of Pa/(K:K).
	double pressureTT;

	/// The second-order partial derivative of pressure with respect to temperature and density in units of Pa/(K:kg/m3).
	double pressureTD;

	/// The second-order partial derivative of pressure with respect to density in units of Pa/((kg/m3):(kg/m3)).
	double pressureDD;
};

/**
 * Outputs the thermodynamic state of water
 */
auto operator<<(std::ostream& out, const WaterThermo& ts) -> std::ostream&;

/**
 * Calculates the thermodynamic state of water
 *
 * It uses the thermodynamic model of Wagner and Pruss (1995) for the calculation.
 *
 * @param T The temperature of water (in units of K)
 * @param P The pressure of water (in units of Pa)
 *
 * @return The thermodynamic state of water, as an instance of @ref WaterThermo
 */
auto waterThermo(double T, double P) -> WaterThermo;

/**
 * Calculates the thermodynamic state of water
 *
 * @param T The temperature of water (in units of K)
 * @param P The pressure of water (in units of Pa)
 * @param model The thermodynamic model as an instance of @ref WaterThermoModel
 *
 * @return The thermodynamic state of water, as an instance of @ref WaterThermo
 */
auto waterThermo(double T, double P, WaterThermoModel model) -> WaterThermo;

/**
 * Calculates the thermodynamic state of water
 *
 * This is a general method that uses the Helmholtz free enery state
 * of water, as an instance of @ref WaterHelmholtz, to completely
 * resolve the thermodynamic state of water.
 *
 * @param T The temperature of water (in units of K)
 * @param P The pressure of water (in units of Pa)
 * @param D The density of water (in units of kg/m3)
 * @param wh The Helmholtz free energy state of water, as an instance of @ref WaterHelmholtz
 *
 * @return The thermodynamic state of water, as an instance of @ref WaterThermo
 */
auto waterThermo(double T, double P, double D, const WaterHelmholtz& wh) -> WaterThermo;

/**
 * Calculates the saturated pressure of water (in units of Pa)
 *
 * It uses the thermodynamic model of Wagner and Pruss (1995) for the calculation.
 *
 * @param T The temperature of water (in units of K)
 *
 * @return The saturated pressure of water (in units of Pa)
 */
auto saturatedPressureWater(double T) -> double;

/**
 * Calculates the saturated liquid-density of water (in units of kg/m3)
 *
 * It uses the thermodynamic model of Wagner and Pruss (1995) for the calculation.
 *
 * @param T The temperature of water (in units of K)
 *
 * @return The saturated liquid-density of water (in units of kg/m3)
 */
auto saturatedLiquidDensityWater(double T) -> double;

/**
 * Calculates the saturated vapour-density of water (in units of kg/m3)
 *
 * It uses the thermodynamic model of Wagner and Pruss (1995) for the calculation.
 *
 * @param T The temperature of water (in units of K)
 *
 * @return The saturated vapour-density of water (in units of kg/m3)
 */
auto saturatedVapourDensityWater(double T) -> double;

/**
 * Calculates the density of water (in units of kg/m3)
 *
 * It uses the thermodynamic model of Wagner and Pruss (1995) for the calculation.
 *
 * @param T The temperature of water (in units of K)
 * @param D The pressure of water (in units of Pa)
 *
 * @return The density of water (in units of kg/m3)
 */
auto densityWater(double T, double P) -> double;

/**
 * Calculates the density of water (in units of kg/m3)
 *
 * @param T The temperature of water (in units of K)
 * @param D The pressure of water (in units of Pa)
 * @param model The thermodynamic model as an instance of @ref WaterThermoModel
 *
 * @return The density of water (in units of kg/m3)
 */
auto densityWater(double T, double P, WaterThermoModel model) -> double;

/**
 * Calculates the pressure of water (in units of Pa)
 *
 * It uses the thermodynamic model of Wagner and Pruss (1995) for the calculation.
 *
 * @param T The temperature of water (in units of K)
 * @param D The density of water (in units of kg/m3)
 *
 * @return The pressure of water (in units of Pa)
 */
auto pressureWater(double T, double D) -> double;

/**
 * Calculates the pressure of water (in units of Pa)
 *
 * @param T The temperature of water (in units of K)
 * @param D The density of water (in units of kg/m3)
 * @param model The thermodynamic model as an instance of @ref WaterThermoModel
 *
 * @return The pressure of water (in units of Pa)
 */
auto pressureWater(double T, double D, WaterThermoModel model) -> double;

} // namespace Reaktor
