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

#include "WaterThermoStateUtils.hpp"

// C++ includes
#include <cmath>
#include <sstream>

// Reaktor includes
#include <Reaktor/Common/Exception.hpp>
#include <Reaktor/Thermodynamics/WaterConstants.hpp>
#include <Reaktor/Thermodynamics/WaterHelmholtzState.hpp>
#include <Reaktor/Thermodynamics/WaterHelmholtzStateHGK.hpp>
#include <Reaktor/Thermodynamics/WaterHelmholtzStateWagnerPruss.hpp>
#include <Reaktor/Thermodynamics/WaterThermoState.hpp>
#include <Reaktor/Thermodynamics/WaterUtils.hpp>

namespace Reaktor {

auto waterThermoStateHGK(double T, double P) -> WaterThermoState
{
    const double D = waterDensityHGK(T, P);
    const WaterHelmholtzState whs = waterHelmholtzStateHGK(T, D);
    return waterThermoState(T, P, D, whs);
}

auto waterThermoStateWagnerPruss(double T, double P) -> WaterThermoState
{
    const double D = waterDensityWagnerPruss(T, P);
    const WaterHelmholtzState whs = waterHelmholtzStateWagnerPruss(T, D);
    return waterThermoState(T, P, D, whs);
}

auto waterThermoState(double T, double P, double D, const WaterHelmholtzState& wh) -> WaterThermoState
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


} // namespace Reaktor
