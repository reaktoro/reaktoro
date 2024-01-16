// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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
#include <Reaktoro/Water/WaterHelmholtzProps.hpp>
#include <Reaktoro/Water/WaterHelmholtzPropsHGK.hpp>
#include <Reaktoro/Water/WaterHelmholtzPropsWagnerPruss.hpp>
#include <Reaktoro/Water/WaterInterpolation.hpp>
#include <Reaktoro/Water/WaterThermoProps.hpp>
#include <Reaktoro/Water/WaterUtils.hpp>

namespace Reaktoro {
namespace {

/// Return a memoized function that computes thermodynamic properties of water using HGK (1984) model.
auto createMemoizedWaterThermoPropsFnHGK()
{
    Fn<WaterThermoProps(real const&, real const&, StateOfMatter)> fn = [](real const& T, real const& P, StateOfMatter som)
    {
        return waterThermoPropsHGK(T, P, som);
    };
    return memoizeLast(fn);
}

/// Return a memoized function that computes thermodynamic properties of water using Wagner & Pruss (1999) model.
auto createMemoizedWaterThermoPropsFnWagnerPruss()
{
    Fn<WaterThermoProps(real const&, real const&, StateOfMatter)> fn = [](real const& T, real const& P, StateOfMatter som)
    {
        return waterThermoPropsWagnerPruss(T, P, som);
    };
    return memoizeLast(fn);
}

/// Return a memoized function that computes thermodynamic properties of water using interpolation on Wagner & Pruss (1999) pre-computed properties.
auto createMemoizedWaterThermoPropsFnWagnerPrussInterp()
{
    Fn<WaterThermoProps(real const&, real const&, StateOfMatter)> fn = [](real const& T, real const& P, StateOfMatter som)
    {
        return waterThermoPropsWagnerPrussInterp(T, P, som);
    };
    return memoizeLast(fn);
}

} // namespace

auto waterThermoPropsHGK(real const& T, real const& P, StateOfMatter som) -> WaterThermoProps
{
    const real D = waterDensityHGK(T, P, som);
    const WaterHelmholtzProps whp = waterHelmholtzPropsHGK(T, D);
    return waterThermoProps(T, P, whp);
}

auto waterThermoPropsHGKMemoized(real const& T, real const& P, StateOfMatter som) -> WaterThermoProps
{
    static thread_local auto fn = createMemoizedWaterThermoPropsFnHGK();
    return fn(T, P, som);
}

auto waterThermoPropsWagnerPruss(real const& T, real const& P, StateOfMatter som) -> WaterThermoProps
{
    const real D = waterDensityWagnerPruss(T, P, som);
    const WaterHelmholtzProps whp = waterHelmholtzPropsWagnerPruss(T, D);
    return waterThermoProps(T, P, whp);
}

auto waterThermoPropsWagnerPrussMemoized(real const& T, real const& P, StateOfMatter som) -> WaterThermoProps
{
    static thread_local auto fn = createMemoizedWaterThermoPropsFnWagnerPruss();
    return fn(T, P, som);
}

auto waterThermoPropsWagnerPrussInterpMemoized(real const& T, real const& P, StateOfMatter som) -> WaterThermoProps
{
    static thread_local auto fn = createMemoizedWaterThermoPropsFnWagnerPrussInterp();
    return fn(T, P, som);
}

auto waterThermoProps(real const& T, real const& P, WaterHelmholtzProps const& whp) -> WaterThermoProps
{
    WaterThermoProps wt;

    // Calculate water density using relation P = \rho^{2}\left(\frac{\partial f}{\partial\rho}\right)_{T}
    auto D = sqrt(P/whp.helmholtzD);

    // Set the temperature of the thermodynamic state of water
    wt.T = T;

    // Set the pressure and its partial derivatives of the thermodynamic state of water
    // wt.P   = P;
    wt.P   = D*D*whp.helmholtzD;
    wt.PD  = 2*D*whp.helmholtzD + D*D*whp.helmholtzDD;
    wt.PT  = D*D*whp.helmholtzTD;
    wt.PDD = 2*whp.helmholtzD + 4*D*whp.helmholtzDD + D*D*whp.helmholtzDDD;
    wt.PTD = 2*D*whp.helmholtzTD + D*D*whp.helmholtzTDD;
    wt.PTT = D*D*whp.helmholtzTTD;

    // Set the density and its partial derivatives of the thermodynamic state of water
    wt.D   = D;
    wt.DT  = -wt.PT/wt.PD;
    wt.DP  =  1.0/wt.PD;
    wt.DTT = -wt.DT*wt.DP*(wt.DT*wt.PDD + 2*wt.PTD + wt.PTT/wt.DT);
    wt.DTP = -wt.DP*wt.DP*(wt.DT*wt.PDD + wt.PTD);
    wt.DPP = -wt.DP*wt.DP*wt.DP*wt.PDD;

    // Set the specific volume of water
    wt.V = 1/D;

    // Set the specific entropy of water
    wt.S = -whp.helmholtzT;

    // Set the specific Helmholtz free energy of water
    wt.A = whp.helmholtz;

    // Set the specific internal energy of water
    wt.U = wt.A + T * wt.S;

    // Set the specific enthalpy of water
    wt.H = wt.U + P/D;

    // Set the specific Gibbs free energy of water
    wt.G = wt.H - T * wt.S;

    // Set the specific isochoric heat capacity of water
    wt.Cv = -T * whp.helmholtzTT;

    // Set the specific isobaric heat capacity of water
    wt.Cp = wt.Cv + T/(D*D)*wt.PT*wt.PT/wt.PD;

    return wt;
}

} // namespace Reaktoro
