// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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

#include "WaterThermoPropsUtils.hpp"

// C++ includes
#include <cmath>
using std::sqrt;

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Memoization.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterHelmholtzProps.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterHelmholtzPropsHGK.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterHelmholtzPropsWagnerPruss.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterThermoProps.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterUtils.hpp>

namespace Reaktoro {
namespace {

/// Return a memoized function that computes thermodynamic properties of water using HGK (1984) model.
auto createMemoizedWaterThermoPropsFnHGK()
{
    Fn<WaterThermoProps(const real&, const real&, StateOfMatter)> fn = [](const real& T, const real& P, StateOfMatter som)
    {
        return waterThermoPropsHGK(T, P, som);
    };
    return memoizeLast(fn);
}

/// Return a memoized function that computes thermodynamic properties of water using Wagner & Pruss (1999) model.
auto createMemoizedWaterThermoPropsFnWagnerPruss()
{
    Fn<WaterThermoProps(const real&, const real&, StateOfMatter)> fn = [](const real& T, const real& P, StateOfMatter som)
    {
        return waterThermoPropsWagnerPruss(T, P, som);
    };
    return memoizeLast(fn);
}

} // namespace

auto waterThermoPropsHGK(real T, real P, StateOfMatter stateofmatter) -> WaterThermoProps
{
    const real D = waterDensityHGK(T, P, stateofmatter);
    const WaterHelmholtzProps whp = waterHelmholtzPropsHGK(T, D);
    return waterThermoProps(T, P, whp);
}

auto waterThermoPropsHGKMemoized(real T, real P, StateOfMatter stateofmatter) -> WaterThermoProps
{
    static thread_local auto fn = createMemoizedWaterThermoPropsFnHGK();
    return fn(T, P, stateofmatter);
}

auto waterThermoPropsWagnerPruss(real T, real P, StateOfMatter stateofmatter) -> WaterThermoProps
{
    const real D = waterDensityWagnerPruss(T, P, stateofmatter);
    const WaterHelmholtzProps whp = waterHelmholtzPropsWagnerPruss(T, D);
    return waterThermoProps(T, P, whp);
}

auto waterThermoPropsWagnerPrussMemoized(real T, real P, StateOfMatter stateofmatter) -> WaterThermoProps
{
    static thread_local auto fn = createMemoizedWaterThermoPropsFnWagnerPruss();
    return fn(T, P, stateofmatter);
}

auto waterThermoProps(real T, real P, const WaterHelmholtzProps& whp) -> WaterThermoProps
{
    WaterThermoProps wt;

    // Calculate water density using relation P = \rho^{2}\left(\frac{\partial f}{\partial\rho}\right)_{T}
    auto D = sqrt(P/whp.helmholtzD);

    // Set the temperature of the thermodynamic state of water
    wt.temperature = T;

    // Set the pressure and its partial derivatives of the thermodynamic state of water
    // wt.pressure   = P;
    wt.pressure   = D*D*whp.helmholtzD;
    wt.pressureD  = 2*D*whp.helmholtzD + D*D*whp.helmholtzDD;
    wt.pressureT  = D*D*whp.helmholtzTD;
    wt.pressureDD = 2*whp.helmholtzD + 4*D*whp.helmholtzDD + D*D*whp.helmholtzDDD;
    wt.pressureTD = 2*D*whp.helmholtzTD + D*D*whp.helmholtzTDD;
    wt.pressureTT = D*D*whp.helmholtzTTD;

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
    wt.entropy = -whp.helmholtzT;

    // Set the specific Helmholtz free energy of water
    wt.helmholtz = whp.helmholtz;

    // Set the specific internal energy of water
    wt.internal_energy = wt.helmholtz + T * wt.entropy;

    // Set the specific enthalpy of water
    wt.enthalpy = wt.internal_energy + P/D;

    // Set the specific Gibbs free energy of water
    wt.gibbs = wt.enthalpy - T * wt.entropy;

    // Set the specific isochoric heat capacity of water
    wt.cv = -T * whp.helmholtzTT;

    // Set the specific isobaric heat capacity of water
    wt.cp = wt.cv + T/(D*D)*wt.pressureT*wt.pressureT/wt.pressureD;

    return wt;
}

} // namespace Reaktoro
