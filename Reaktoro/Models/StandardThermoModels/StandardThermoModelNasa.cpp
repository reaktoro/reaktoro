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

#include "StandardThermoModelNasa.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Serialization/Models/StandardThermoModels.hpp>

namespace Reaktoro {
namespace detail {

auto indexTemperatureInterval(const Vec<StandardThermoModelParamsNasa::Polynomial>& polynomials, const real& T) -> Index
{
    for(auto i = 0; i < polynomials.size(); ++i)
    {
        const auto Tmin = polynomials[i].Tmin;
        const auto Tmax = polynomials[i].Tmax;
        if(Tmin <= T && T <= Tmax)
            return i;
    }
    return polynomials.size();
}

auto computeStandardThermoProps(const StandardThermoModelParamsNasa::Polynomial& polynomial, const real& T) -> StandardThermoProps
{
    StandardThermoProps props;

    const auto& Tmin = polynomial.Tmin;
    const auto& Tmax = polynomial.Tmax;

    assert(Tmin <= T && T <= Tmax);

    const auto& a1 = polynomial.a1;
    const auto& a2 = polynomial.a2;
    const auto& a3 = polynomial.a3;
    const auto& a4 = polynomial.a4;
    const auto& a5 = polynomial.a5;
    const auto& a6 = polynomial.a6;
    const auto& a7 = polynomial.a7;
    const auto& b1 = polynomial.b1;
    const auto& b2 = polynomial.b2;

    const auto T2 = T*T;
    const auto T3 = T*T2;
    const auto T4 = T*T3;
    const auto lnT = log(T);

    const auto R = universalGasConstant;

    const auto Cp0 = ( a1/T2 + a2/T + a3 + a4*T + a5*T2 + a6*T3 + a7*T4) * R;
    const auto H0  = (-a1/T2 + a2*lnT/T + a3 + a4*T/2.0 + a5*T2/3.0 + a6*T3/4.0 + a7*T4/5.0 + b1/T) * R*T;
    const auto S0  = (-a1/T2*0.5 - a2/T + a3*lnT + a4*T + a5*T2/2.0 + a6*T3/3.0 + a7*T4/4.0 + b2) * R;

    props.G0  = H0 - T*S0;
    props.H0  = H0;
    props.V0  = 0.0;
    props.VT0 = 0.0;
    props.VP0 = 0.0;
    props.Cp0 = Cp0;

    return props;
}

auto computeStandardThermoProps(const StandardThermoModelParamsNasa& params, const real& T) -> StandardThermoProps
{
    // Check if params corresponds to a species without temperature intervals, and just enthalpy at a single temperature point.
    if(params.polynomials.empty())
    {
        StandardThermoProps props;
        props.G0 = params.H0; // NOTE: No given data for computation of G0, so assuming G0 = H0 at T0
        props.H0 = params.H0;
        return props;
    };

    // Find the index of the temperature interval in which T is contained
    const auto iT = indexTemperatureInterval(params.polynomials, T);

    // Compute the standard thermodynamic properties only if within valid temperature range
    if(iT < params.polynomials.size())
        return computeStandardThermoProps(params.polynomials[iT], T);

    // Otherwise, return standard thermo props whose G0 is high to penalize the species from appearing at equilibrium
    StandardThermoProps props;
    props.G0 = 999'999'999'999; // NOTE: This high value for G0 is to ensure the condensed species with limited data is never stable at equilibrium
    return props;
}

} // namemespace detail

auto StandardThermoModelNasa(const StandardThermoModelParamsNasa& params) -> StandardThermoModel
{
    auto evalfn = [=](StandardThermoProps& props, real T, real P)
    {
        props = detail::computeStandardThermoProps(params, T);
    };

    Data paramsdata;
    paramsdata["Nasa"] = params;

    return StandardThermoModel(evalfn, paramsdata);
}

} // namespace Reaktoro
