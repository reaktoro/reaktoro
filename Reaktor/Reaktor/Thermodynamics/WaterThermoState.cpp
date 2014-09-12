// Reaktor is a C++ library for computational reaction modelling.
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

#include "WaterThermoState.hpp"

// C++ includes
#include <cmath>
#include <sstream>

// Reaktor includes
#include <Reaktor/Common/Exception.hpp>
#include <Reaktor/Thermodynamics/WaterConstants.hpp>
#include <Reaktor/Thermodynamics/WaterHelmholtzState.hpp>
#include <Reaktor/Thermodynamics/WaterThermoStateWagnerPruss.hpp>

namespace Reaktor {

WaterThermoState::WaterThermoState()
: volume(0), entropy(0), helmholtz(0), internal_energy(0), enthalpy(0), gibbs(0), cv(0), cp(0),
  density(0), densityT(0), densityP(0), densityTT(0), densityTP(0), densityPP(0),
  pressureT(0), pressureD(0), pressureTT(0), pressureTD(0), pressureDD(0)
{}

auto operator<<(std::ostream& out, const WaterThermoState& ts) -> std::ostream&
{
	out << "volume     = " << ts.volume     << std::endl;
	out << "entropy    = " << ts.entropy    << std::endl;
	out << "helmholtz  = " << ts.helmholtz  << std::endl;
	out << "energy     = " << ts.internal_energy     << std::endl;
	out << "enthalpy   = " << ts.enthalpy   << std::endl;
	out << "gibbs      = " << ts.gibbs      << std::endl;
	out << "cv         = " << ts.cv         << std::endl;
	out << "cp         = " << ts.cp         << std::endl;
	out << "density    = " << ts.density    << std::endl;
	out << "densityT   = " << ts.densityT   << std::endl;
	out << "densityP   = " << ts.densityP   << std::endl;
	out << "densityTT  = " << ts.densityTT  << std::endl;
	out << "densityTP  = " << ts.densityTP  << std::endl;
	out << "densityPP  = " << ts.densityPP  << std::endl;
	out << "pressureT  = " << ts.pressureT  << std::endl;
	out << "pressureD  = " << ts.pressureD  << std::endl;
	out << "pressureTT = " << ts.pressureTT << std::endl;
	out << "pressureTD = " << ts.pressureTD << std::endl;
	out << "pressureDD = " << ts.pressureDD << std::endl;

	return out;
}

auto waterThermo(double T, double P) -> WaterThermoState
{
	return waterThermo(T, P, WagnerPruss);
}

auto waterThermo(double T, double P, WaterThermoModel model) -> WaterThermoState
{
	const double D = densityWater(T, P, model);

	WaterHelmholtzState wh = waterHelmholtz(T, D, model);

	return waterThermo(T, P, D, wh);
}

auto waterThermo(double T, double P, double D, const WaterHelmholtzState& wh) -> WaterThermoState
{
	WaterThermoState wt;

	// Set the temperature of the thermodynamic state of water
	wt.temperature = T;

	// Set the pressure and its partial derivatives of the thermodynamic state of water
	wt.pressure   = P;
	wt.pressureD  = 2*D*wh.helmholtzD + D*D*wh.helmholtzDD;
	wt.pressureT  = D*D*wh.helmholtzTD;
	wt.pressureDD = 2*wh.helmholtzD + 4*D*wh.helmholtzDD + D*D*wh.helmholtzDDD;
	wt.pressureTD = 2*D*wh.helmholtzTD + D*D*wh.helmholtzTDD;
	wt.pressureTT = D*D*wh.helmholtzTTD;

	// Set the density and its partial derivatives of the thermodynamic state of water
	wt.density   = D;
	wt.densityT  = -wt.pressureT/wt.pressureD;
	wt.densityP  =  1.0/wt.pressureD;
	wt.densityTT = -wt.densityT*wt.densityP*(wt.densityT*wt.pressureDD + 2*wt.pressureTD + wt.pressureTT/wt.densityT);
	wt.densityTP = -wt.densityP*wt.densityP*(wt.densityT*wt.pressureDD + wt.pressureTD);
	wt.densityPP = -wt.densityP*wt.densityP*wt.densityP*wt.pressureDD;

	// Set the specific volume of water
	wt.volume = 1/D;

	// Set the specific entropy of water
	wt.entropy = -wh.helmholtzT;

	// Set the specific Helmholtz free energy of water
	wt.helmholtz = wh.helmholtz;

	// Set the specific internal energy of water
	wt.internal_energy = wt.helmholtz + T * wt.entropy;

	// Set the specific enthalpy of water
	wt.enthalpy = wt.internal_energy + P/D;

	// Set the specific Gibbs free energy of water
	wt.gibbs = wt.enthalpy - T * wt.entropy;

	// Set the specific isochoric heat capacity of water
	wt.cv = -T * wh.helmholtzTT;

	// Set the specific isobaric heat capacity of water
	wt.cp = wt.cv + T/(D*D)*wt.pressureT*wt.pressureT/wt.pressureD;

	return wt;
}

auto saturatedPressureWater(double T) -> double
{
	return saturatedPressureWaterWagnerPruss(T);
}

auto saturatedLiquidDensityWater(double T) -> double
{
	return saturatedLiquidDensityWaterWagnerPruss(T);
}

auto saturatedVapourDensityWater(double T) -> double
{
	return saturatedVapourDensityWaterWagnerPruss(T);
}

auto densityWater(double T, double P) -> double
{
	return densityWater(T, P, WagnerPruss);
}

auto densityWater(double T, double P, WaterThermoModel model) -> double
{
	// Auxiliary constants for the Newton's iterations
	const int max_iters = 100;
	const double tolerance = 1.0e-08;

	// Determine the physical state of water, where: 0-liquid, 1-vapour, 2-supercritical
	int state;

	if(T <= waterCriticalTemperature)
		state = (P <= saturatedPressureWater(T)) ? 1 : 0;
	else
		state = 2;

	// Determine an adequate initial guess for (dimensionless) density based on the physical state of water
	double D; switch(state)
	{
	case 0: D = saturatedLiquidDensityWater(T); break;
	case 1: D = saturatedVapourDensityWater(T); break;
	case 2: D = waterCriticalDensity * 0.99; break; // some derivatives of the Helmholtz free energy do not exist in the critical point
	}

	// Apply the Newton's method to the pressure-density equation
	for(int i = 1; i <= max_iters; ++i)
	{
		WaterHelmholtzState h = waterHelmholtz(T, D, model);

		const double f  = (D*D*h.helmholtzD - P)/waterCriticalPressure;
		const double df = (2*D*h.helmholtzD + D*D*h.helmholtzDD)/waterCriticalPressure;

		D = (D > f/df) ? D - f/df : P/(D*h.helmholtzD);

		if(std::abs(f) < tolerance)
			return D;
	}

	Exception exception;
    exception.error << "Unable to calculate the density of water.";
    exception.reason << "The calculations did not converge at temperature " << T << " K and pressure " << P << "Pa.";
    raise(exception);

	return 0;
}

auto pressureWater(double T, double D) -> double
{
	return pressureWater(T, D, WagnerPruss);
}

auto pressureWater(double T, double D, WaterThermoModel model) -> double
{
	WaterHelmholtzState h = waterHelmholtz(T, D, model);

	return D*D*h.helmholtzD;
}

} // namespace Reaktor
