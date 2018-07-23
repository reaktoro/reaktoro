// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/ThermoScalar.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterHelmholtzState.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterHelmholtzStateHGK.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterHelmholtzStateWagnerPruss.hpp>

namespace Reaktoro {

template<typename HelmholtsModel>
auto waterDensity(Temperature T, Pressure P, const HelmholtsModel& model, StateOfMatter stateofmatter) -> ThermoScalar
{
    // Auxiliary constants for the Newton's iterations
    const auto max_iters = 100;
    const auto tolerance = 1.0e-08;

    // Auxiliary constants for initial guess computation for water density
    const auto R = universalGasConstant;
    const auto Twc = waterCriticalTemperature;
    const auto Pwc = waterCriticalPressure;

    // Compute initial guess for molar volume (m3/mol) according to suggestion in Smith,
    const auto Vwc_liquid = R*Twc/Pwc; // initial guess for volume if liquid state (see equation 3.55, page 97 in Smith, Van Ness, Abbott, 2005)
    const auto Vwc_vapor = R*T/P; // initial guess for volume if gas state (see equation 3.49, page 96 in Smith, Van Ness, Abbott, 2005)

    // Compute initial guess for density (kg/m3)
//    const auto Dwc_liquid = waterMolarMass / Vwc_liquid;
    const auto Dwc_liquid = 10.0 * waterCriticalDensity;
    const auto Dwc_vapor = waterMolarMass / Vwc_vapor;

    // Determine an adequate initial guess for (dimensionless) density based on the physical state of water
    ThermoScalar D;

    switch(stateofmatter)
    {
    case StateOfMatter::Liquid: D = Dwc_liquid; break;
    default: D = Dwc_vapor; break;
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

auto waterDensityHGK(Temperature T, Pressure P, StateOfMatter stateofmatter) -> ThermoScalar
{
    return waterDensity(T, P, waterHelmholtzStateHGK, stateofmatter);
}

auto waterLiquidDensityHGK(Temperature T, Pressure P) -> ThermoScalar
{
    return waterDensityHGK(T, P, StateOfMatter::Liquid);
}

auto waterVaporDensityHGK(Temperature T, Pressure P) -> ThermoScalar
{
    return waterDensityHGK(T, P, StateOfMatter::Gas);
}

auto waterDensityWagnerPruss(Temperature T, Pressure P, StateOfMatter stateofmatter) -> ThermoScalar
{
    return waterDensity(T, P, waterHelmholtzStateWagnerPruss, stateofmatter);
}

auto waterLiquidDensityWagnerPruss(Temperature T, Pressure P) -> ThermoScalar
{
    return waterDensityWagnerPruss(T, P, StateOfMatter::Liquid);
}

auto waterVaporDensityWagnerPruss(Temperature T, Pressure P) -> ThermoScalar
{
    return waterDensityWagnerPruss(T, P, StateOfMatter::Gas);
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
