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

#include "WaterUtils.hpp"

// C++ includes
#include <cmath>
using std::abs;
using std::exp;
using std::pow;

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Water/WaterConstants.hpp>
#include <Reaktoro/Water/WaterHelmholtzProps.hpp>
#include <Reaktoro/Water/WaterHelmholtzPropsHGK.hpp>
#include <Reaktoro/Water/WaterHelmholtzPropsWagnerPruss.hpp>
#include <Reaktoro/Water/WaterInterpolation.hpp>

namespace Reaktoro {

template<typename HelmholtsModel>
auto waterDensity(real const& T, real const& P, HelmholtsModel const& model, StateOfMatter stateofmatter) -> real
{
    // Auxiliary constants for the Newton's iterations
    const auto max_iters = 100;
    const auto tolerance = 1.0e-06;

    // Determine an adequate initial guess for density based on the desired physical state of water
    real D = waterDensityWagnerPrussInterp(T, P, stateofmatter);

    for(int i = 1; i <= max_iters; ++i)
    {
        WaterHelmholtzProps h = model(T, D);

        const auto AD = h.helmholtzD;
        const auto ADD = h.helmholtzDD;
        const auto ADDD = h.helmholtzDDD;

        const auto F = D*D*AD/P - 1;
        const auto FD = (2*D*AD + D*D*ADD)/P;
        const auto FDD = (2*AD + 2*D*ADD + 2*D*ADD + D*D*ADDD)/P;

        const auto f = 0.5 * F*F;
        const auto g = F*FD;
        const auto H = FD*FD + F*FDD;

        if(D > g/H)
            D -= g/H;
        else if(D > F/FD)
            D -= F/FD;
        else D *= 0.1;

        if(abs(F) < tolerance || abs(g) < tolerance)
            return D;
    }

    errorif(true, "Unable to calculate the density of water because the calculations did not converge at temperature ", T, " K and pressure ", P, " Pa.");

    return {};
}

auto waterDensityHGK(real const& T, real const& P, StateOfMatter stateofmatter) -> real
{
    return waterDensity(T, P, waterHelmholtzPropsHGK, stateofmatter);
}

auto waterLiquidDensityHGK(real const& T, real const& P) -> real
{
    return waterDensityHGK(T, P, StateOfMatter::Liquid);
}

auto waterVaporDensityHGK(real const& T, real const& P) -> real
{
    return waterDensityHGK(T, P, StateOfMatter::Gas);
}

auto waterDensityWagnerPruss(real const& T, real const& P, StateOfMatter stateofmatter) -> real
{
    return waterDensity(T, P, waterHelmholtzPropsWagnerPruss, stateofmatter);
}

auto waterLiquidDensityWagnerPruss(real const& T, real const& P) -> real
{
    return waterDensityWagnerPruss(T, P, StateOfMatter::Liquid);
}

auto waterVaporDensityWagnerPruss(real const& T, real const& P) -> real
{
    return waterDensityWagnerPruss(T, P, StateOfMatter::Gas);
}

template<typename HelmholtzModel>
auto waterPressure(real const& T, real const& D, HelmholtzModel const& model) -> real
{
    WaterHelmholtzProps h = model(T, D);
    return D*D*h.helmholtzD;
}

auto waterPressureHGK(real const& T, real const& D) -> real
{
    return waterPressure(T, D, waterHelmholtzPropsHGK);
}

auto waterPressureWagnerPruss(real const& T, real const& D) -> real
{
    return waterPressure(T, D, waterHelmholtzPropsHGK);
}

auto waterSaturationPressureWagnerPruss(real const& T) -> real
{
    const double a1 = -7.85951783;
    const double a2 =  1.84408259;
    const double a3 = -11.7866497;
    const double a4 =  22.6807411;
    const double a5 = -15.9618719;
    const double a6 =  1.80122502;

    const double Tcr = waterCriticalTemperature;
    const double Pcr = waterCriticalPressure;

    const auto t   = 1 - T/Tcr;
    const auto t15 = pow(t, 1.5);
    const auto t30 = t15 * t15;
    const auto t35 = t15 * t * t;
    const auto t40 = t30 * t;
    const auto t75 = t35 * t40;

    return Pcr * exp(Tcr/T * (a1*t + a2*t15 + a3*t30 + a4*t35 + a5*t40 + a6*t75));
}

auto waterSaturationLiquidDensityWagnerPruss(real const& T) -> real
{
    const double b1 =  1.99274064;
    const double b2 =  1.09965342;
    const double b3 = -0.510839303;
    const double b4 = -1.75493479;
    const double b5 = -45.5170352;
    const double b6 = -6.74694450e+05;

    const double Tcr = waterCriticalTemperature;
    const double Dcr = waterCriticalDensity;

    const auto t     = 1 - T/Tcr;
    const auto t13   = pow(t, 1./3);
    const auto t23   = t13 * t13;
    const auto t53   = t13 * t23 * t23;
    const auto t163  = t13 * t53 * t53 * t53;
    const auto t433  = t163 * t163 * t53 * t * t;
    const auto t1103 = t433 * t433 * t163 * t53 * t;

    return Dcr * (1 + b1*t13 + b2*t23 + b3*t53 + b4*t163 + b5*t433 + b6*t1103);
}

auto waterSaturationVapourDensityWagnerPruss(real const& T) -> real
{
    const double c1 = -2.03150240;
    const double c2 = -2.68302940;
    const double c3 = -5.38626492;
    const double c4 = -17.2991605;
    const double c5 = -44.7586581;
    const double c6 = -63.9201063;

    const double Tcr = waterCriticalTemperature;
    const double Dcr = waterCriticalDensity;

    const auto t    = 1 - T/Tcr;
    const auto t16  = pow(t, 1./6);
    const auto t26  = t16 * t16;
    const auto t46  = t26 * t26;
    const auto t86  = t46 * t46;
    const auto t186 = t86 * t86 * t26;
    const auto t376 = t186 * t186 * t16;
    const auto t716 = t376 * t186 * t86 * t86;

    return Dcr * exp(Tcr/T * (c1*t26 + c2*t46 + c3*t86 + c4*t186 + c5*t376 + c6*t716));
}

auto waterSaturatedPressureWagnerPruss(real const& T) -> real
{
    errorif(true, "waterSaturatedPressureWagnerPruss has been renamed to waterSaturationPressureWagnerPruss")
}

auto waterSaturatedLiquidDensityWagnerPruss(real const& T) -> real
{
    errorif(true, "waterSaturatedLiquidDensityWagnerPruss has been renamed to waterSaturationLiquidDensityWagnerPruss")
}

auto waterSaturatedVapourDensityWagnerPruss(real const& T) -> real
{
    errorif(true, "waterSaturatedVapourDensityWagnerPruss has been renamed to waterSaturationVapourDensityWagnerPruss")
}

} // namespace Reaktoro
