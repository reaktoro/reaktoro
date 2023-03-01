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

#include "StandardThermoModelHKF.hpp"

// C++ includes
#include <cmath>
using std::log;

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Memoization.hpp>
#include <Reaktoro/Models/StandardThermoModels/Support/SpeciesElectroProps.hpp>
#include <Reaktoro/Models/StandardThermoModels/Support/SpeciesElectroPropsHKF.hpp>
#include <Reaktoro/Serialization/Models/StandardThermoModels.hpp>
#include <Reaktoro/Water/WaterElectroProps.hpp>
#include <Reaktoro/Water/WaterElectroPropsJohnsonNorton.hpp>
#include <Reaktoro/Water/WaterInterpolation.hpp>
#include <Reaktoro/Water/WaterThermoProps.hpp>
#include <Reaktoro/Water/WaterThermoPropsUtils.hpp>

namespace Reaktoro {
namespace {

/// The reference temperature assumed in the HKF equations of state (in units of K)
const auto Tr = 298.15;

/// The reference pressure assumed in the HKF equations of state (in units of Pa)
const auto Pr = 1.0e+05;

/// The reference Born function Z (dimensionless)
const auto Zr = -1.278055636e-02;

/// The reference Born function Y (dimensionless)
const auto Yr = -5.795424563e-05;

/// The constant characteristics @eq{\Theta} of the solvent (in units of K)
const auto theta = 228.0;

/// The constant characteristics @eq{\Psi} of the solvent (in units of Pa)
const auto psi = 2600.0e+05;

/// Return a memoized function that computes thermodynamic properties of water using Wagner & Pruss (1999) model.
auto createMemoizedWaterElectroPropsFnJohnsonNorton()
{
    Fn<WaterElectroProps(const real&, const real&)> fn = [](const real& T, const real& P)
    {
        const auto wtp = waterThermoPropsWagnerPrussMemoized(T, P, StateOfMatter::Liquid);
        return Reaktoro::waterElectroPropsJohnsonNorton(T, P, wtp);
    };
    return memoizeLast(fn);
}

/// Return the computed electrostatic properties of water at @p T and @p P using Johnson and Norton (1991) model.
auto memoizedWaterElectroPropsJohnsonNorton(const real& T, const real& P) -> WaterElectroProps
{
    static thread_local auto fn = createMemoizedWaterElectroPropsFnJohnsonNorton();
    return fn(T, P);
}

} // namespace

/// Return a Vec<Param> object containing all Param objects in @p params.
auto extractParams(const StandardThermoModelParamsHKF& params) -> Vec<Param>
{
    const auto& [Gf, Hf, Sr, a1, a2, a3, a4, c1, c2, wref, charge, Tmax] = params;
    return {Gf, Hf, Sr, a1, a2, a3, a4, c1, c2, wref};
}

/// Return a ModelSerializer for given model parameters in @p params.
auto createModelSerializer(const StandardThermoModelParamsHKF& params) -> ModelSerializer
{
    return [=]()
    {
        Data node;
        node["HKF"] = params;
        return node;
    };
}

auto StandardThermoModelHKF(const StandardThermoModelParamsHKF& params) -> StandardThermoModel
{
    waterThermoPropsWagnerPrussInterpData(StateOfMatter::Liquid); // this call exists to force an initialization operation so that when waterThermoPropsWagnerPrussInterp is called for the first time, this initialization has been performed already.

    auto evalfn = [=](StandardThermoProps& props, real T, real P)
    {
        auto& [G0, H0, V0, Cp0, VT0, VP0] = props;
        const auto& [Gf, Hf, Sr, a1, a2, a3, a4, c1, c2, wr, charge, Tmax] = params;

        const auto wtp = waterThermoPropsWagnerPrussMemoized(T, P, StateOfMatter::Liquid);
        const auto wep = memoizedWaterElectroPropsJohnsonNorton(T, P);
        const auto gstate = gHKF::compute(T, P, wtp);
        const auto aep = speciesElectroPropsHKF(gstate, params);

        const auto& w   = aep.w;
        const auto& wT  = aep.wT;
        const auto& wP  = aep.wP;
        const auto& wTP = aep.wTP;
        const auto& wTT = aep.wTT;
        const auto& wPP = aep.wPP;
        const auto& Z   = wep.bornZ;
        const auto& Y   = wep.bornY;
        const auto& Q   = wep.bornQ;
        const auto& U   = wep.bornU;
        const auto& N   = wep.bornN;
        const auto& X   = wep.bornX;
        const auto Tth  = T - theta;
        const auto Tth2 = Tth*Tth;
        const auto Tth3 = Tth*Tth2;

        V0 = a1 + a2/(psi + P) + (a3 + a4/(psi + P))/(T - theta) - w*Q - (Z + 1)*wP;

        VT0 = -(a3 + a4/(psi + P))/((T - theta)*(T - theta)) - wT*Q - w*U - Y*wP - (Z + 1)*wTP;

        VP0 = -a2/((psi + P)*(psi + P)) + (-a4/((psi + P)*(psi + P)))/(T - theta) - wP*Q - w*N - Q*wP - (Z + 1)*wPP;

        G0 = Gf - Sr*(T - Tr) - c1*(T*log(T/Tr) - T + Tr)
            + a1*(P - Pr) + a2*log((psi + P)/(psi + Pr))
            - c2*((1.0/(T - theta) - 1.0/(Tr - theta))*(theta - T)/theta
            - T/(theta*theta)*log(Tr/T * (T - theta)/(Tr - theta)))
            + 1.0/(T - theta)*(a3*(P - Pr) + a4*log((psi + P)/(psi + Pr)))
            - w*(Z + 1) + wr*(Zr + 1) + wr*Yr*(T - Tr);

        H0 = Hf + c1*(T - Tr) - c2*(1.0/(T - theta) - 1.0/(Tr - theta))
            + a1*(P - Pr) + a2*log((psi + P)/(psi + Pr))
            + (2.0*T - theta)/Tth2*(a3*(P - Pr)
            + a4*log((psi + P)/(psi + Pr)))
            - w*(Z + 1) + w*T*Y + T*(Z + 1)*wT + wr*(Zr + 1) - wr*Tr*Yr;

        Cp0 = c1 + c2/Tth2 - 2.0*T/Tth3*(a3*(P - Pr) + a4*log((psi + P)/(psi + Pr))) + w*T*X + 2.0*T*Y*wT + T*(Z + 1.0)*wTT;

        // S0 = Sr + c1*log(T/Tr)
        //     - c2/theta*(1.0/(T - theta)
        //     - 1.0/(Tr - theta) + log(Tr/T * (T - theta)/(Tr - theta))/theta)
        //     + 1.0/Tth2*(a3*(P - Pr)
        //     + a4*log((psi + P)/(psi + Pr)))
        //     + w*Y + (Z + 1)*wT - wr*Yr;
    };

    return StandardThermoModel(evalfn, extractParams(params), createModelSerializer(params));
}

} // namespace Reaktoro
