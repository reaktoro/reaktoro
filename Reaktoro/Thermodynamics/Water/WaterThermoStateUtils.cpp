// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/ThermoScalar.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterHelmholtzState.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterHelmholtzStateHGK.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterHelmholtzStateWagnerPruss.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterThermoState.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterUtils.hpp>

namespace Reaktoro {

auto waterThermoStateHGK(Temperature T, Pressure P, StateOfMatter stateofmatter) -> WaterThermoState
{
    const ThermoScalar D = waterDensityHGK(T, P, stateofmatter);
    const WaterHelmholtzState whs = waterHelmholtzStateHGK(T, D);
    return waterThermoState(T, P, whs);
}

auto waterThermoStateWagnerPruss(Temperature T, Pressure P, StateOfMatter stateofmatter) -> WaterThermoState
{
    const ThermoScalar D = waterDensityWagnerPruss(T, P, stateofmatter);
    const WaterHelmholtzState whs = waterHelmholtzStateWagnerPruss(T, D);
    return waterThermoState(T, P, whs);
}

auto waterThermoState(Temperature T, Pressure P, const WaterHelmholtzState& whs) -> WaterThermoState
{
	WaterThermoState wt;

	// Calculate water density using relation P = \rho^{2}\left(\frac{\partial f}{\partial\rho}\right)_{T}
    auto D = sqrt(P/whs.helmholtzD);

	// Set the temperature of the thermodynamic state of water
	wt.temperature = T;

	// Set the pressure and its partial derivatives of the thermodynamic state of water
	// wt.pressure   = P;
	wt.pressure   = D*D*whs.helmholtzD;
	wt.pressureD  = 2*D*whs.helmholtzD + D*D*whs.helmholtzDD;
	wt.pressureT  = D*D*whs.helmholtzTD;
	wt.pressureDD = 2*whs.helmholtzD + 4*D*whs.helmholtzDD + D*D*whs.helmholtzDDD;
	wt.pressureTD = 2*D*whs.helmholtzTD + D*D*whs.helmholtzTDD;
	wt.pressureTT = D*D*whs.helmholtzTTD;

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
	wt.entropy = -whs.helmholtzT;

	// Set the specific Helmholtz free energy of water
	wt.helmholtz = whs.helmholtz;

	// Set the specific internal energy of water
	wt.internal_energy = wt.helmholtz + T * wt.entropy;

	// Set the specific enthalpy of water
	wt.enthalpy = wt.internal_energy + P/D;

	// Set the specific Gibbs free energy of water
	wt.gibbs = wt.enthalpy - T * wt.entropy;

	// Set the specific isochoric heat capacity of water
	wt.cv = -T * whs.helmholtzTT;

	// Set the specific isobaric heat capacity of water
	wt.cp = wt.cv + T/(D*D)*wt.pressureT*wt.pressureT/wt.pressureD;

	return wt;
}

} // namespace Reaktoro
