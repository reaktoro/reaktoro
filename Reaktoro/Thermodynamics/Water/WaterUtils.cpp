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

#include "WaterUtils.hpp"

// C++ includes
#include <cmath>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/ThermoScalar.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterHelmholtzState.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterHelmholtzStateHGK.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterHelmholtzStateWagnerPruss.hpp>

namespace Reaktoro {

template<typename HelmholtsModel>
auto waterDensity(Temperature T, Pressure P, const HelmholtsModel& model) -> ThermoScalar
{
    // Auxiliary constants for the Newton's iterations
    const int max_iters = 100;
    const double tolerance = 1.0e-08;

    // Determine the physical state of water, where: 0-liquid, 1-vapour, 2-supercritical
    int state;

    if(T.val <= waterCriticalTemperature)
        state = (P <= waterSaturatedPressureWagnerPruss(T)) ? 1 : 0;
    else
        state = 2;

    // Determine an adequate initial guess for (dimensionless) density based on the physical state of water
    ThermoScalar D;

    switch(state)
    {
    case 0: D = waterSaturatedLiquidDensityWagnerPruss(T); break;
    case 1: D = waterSaturatedVapourDensityWagnerPruss(T); break;
    case 2: D = waterCriticalDensity * 0.99; break; // some derivatives of the Helmholtz free energy do not exist in the critical point
    }

    // Apply the Newton's method to the pressure-density equation
    for(int i = 1; i <= max_iters; ++i)
    {
        WaterHelmholtzState h = model(T, D);

        const auto f  = (D*D*h.helmholtzD - P)/waterCriticalPressure;
        const auto df = (2*D*h.helmholtzD + D*D*h.helmholtzDD)/waterCriticalPressure;

        D = (D > f/df) ? D - f/df : P/(D*h.helmholtzD);

        if(abs(f) < tolerance)
            return D;
    }

    Exception exception;
    exception.error << "Unable to calculate the density of water.";
    exception.reason << "The calculations did not converge at temperature "
        << T.val << " K and pressure " << P.val << "Pa.";
    RaiseError(exception);

    return {};
}

auto waterDensityHGK(Temperature T, Pressure P) -> ThermoScalar
{
    return waterDensity(T, P, waterHelmholtzStateHGK);
}

auto waterDensityWagnerPruss(Temperature T, Pressure P) -> ThermoScalar
{
    return waterDensity(T, P, waterHelmholtzStateWagnerPruss);
}

template<typename HelmholtzModel>
auto waterPressure(Temperature T, ThermoScalar D, const HelmholtzModel& model) -> ThermoScalar
{
    WaterHelmholtzState h = model(T, D);
    return D*D*h.helmholtzD.val;
}

auto waterPressureHGK(Temperature T, ThermoScalar D) -> ThermoScalar
{
    return waterPressure(T, D, waterHelmholtzStateHGK);
}

auto waterPressureWagnerPruss(Temperature T, ThermoScalar D) -> ThermoScalar
{
    return waterPressure(T, D, waterHelmholtzStateHGK);
}

auto waterSaturatedPressureWagnerPruss(Temperature T) -> ThermoScalar
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

auto waterSaturatedLiquidDensityWagnerPruss(Temperature T) -> ThermoScalar
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

auto waterSaturatedVapourDensityWagnerPruss(Temperature T) -> ThermoScalar
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

} // namespace Reaktoro
