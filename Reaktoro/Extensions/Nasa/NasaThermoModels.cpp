// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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

#include "NasaThermoModels.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>

namespace Reaktoro {

namespace NasaUtils {

/// Return either the original given temperature or one of its lower or upper bounds to stay within valid range.
auto correctTemperature(const real& Torig, const real& Tmin, const real& Tmax, const String& spname) -> real
{
    const auto T = (Torig < Tmin) ? Tmin : (Torig > Tmax) ? Tmax : Torig;

    warning(T != Torig, "Computing standard thermodynamic properties of "
        "species ", spname, " at ", Torig, " K, which is outside "
        "the valid temperature interval of T(min) = ", Tmin, " K "
        "and T(max) = ", Tmax, " K. Proceeding instead with T = ", T, " K "
        "for the computation of these standard properties.");
}

auto getNasaThermoParamsForGivenTemperature(const Vec<NasaThermoParams>& params, const real& T) -> NasaThermoParams
{
    for(const auto& obj : params)
        if(obj.Tmin <= T && T <= obj.Tmax)
            return obj;
    return {};
}

auto computeStandardThermoProps(const NasaThermoParams& params, const real& T) -> StandardThermoProps
{
    StandardThermoProps props;

    const auto& Tmin = params.Tmin;
    const auto& Tmax = params.Tmax;

    assert(Tmin <= T && T <= Tmax);

    const auto& a1 = params.a1;
    const auto& a2 = params.a2;
    const auto& a3 = params.a3;
    const auto& a4 = params.a4;
    const auto& a5 = params.a5;
    const auto& a6 = params.a6;
    const auto& a7 = params.a7;
    const auto& b1 = params.b1;
    const auto& b2 = params.b2;

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
    props.Cp0 = Cp0;
    props.Cv0 = 0.0;

    return props;
}

auto computeStandardThermoProps(const NasaSpecies& species, const real& Torig) -> StandardThermoProps
{
    const auto T = correctTemperature(Torig, species.Tmin, species.Tmax, species.name);
    const auto params = getNasaThermoParamsForGivenTemperature(species.thermodata, T);
    return computeStandardThermoProps(params, T);
}

} // namespace NasaUtils

auto StandardThermoModelNasa(const NasaSpecies& species) -> StandardThermoModel
{
    if(species.thermodata.empty())
        return [species](real T, real P)
        {
            warning(true, "Cannot compute the standard thermodynamic "
                "properties of species ", species.name, " since it does not have "
                "sufficient NASA CEA parameters (i.e., its number of temperature "
                "intervals is zero). Returning zero for all such properties.");
            return StandardThermoProps{};
        };

    return [species](real T, real P)
    {
        return NasaUtils::computeStandardThermoProps(species, T);
    };
}

} // namespace Reaktoro
