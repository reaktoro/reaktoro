// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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

#include "StandardThermoModelWaterHKF.hpp"

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterThermoState.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterThermoStateUtils.hpp>

namespace Reaktoro {

/// Return a Params object containing all Param objects in @p params.
auto extractParams(const StandardThermoModelParamsWaterHKF& params) -> Params
{
    return {}; // there are no Param members in StandardThermoModelParamsWaterHKF
}

auto StandardThermoModelWaterHKF(const StandardThermoModelParamsWaterHKF& params) -> StandardThermoModel
{
    auto evalfn = [=](StandardThermoProps& props, real T, real P)
    {
        auto& [G0, H0, V0, Cp0, Cv0] = props;
        const auto& [Ttr, Str, Gtr, Htr] = params;

        const auto wtp = waterThermoStateWagnerPrussMemoized(T, P, StateOfMatter::Liquid);

        // Convert from specific properties to molar properties
        const auto Sw = waterMolarMass * wtp.entropy;         // from J/(kg*K) to J/(mol*K)
        const auto Hw = waterMolarMass * wtp.enthalpy;        // from J/kg to J/mol
        const auto Uw = waterMolarMass * wtp.internal_energy; // from J/kg to J/mol

        // See Helgeson and Kirkham (1974), page 1098.
        H0  = Hw + Htr;
        G0  = Hw - T*(Sw + Str) + Ttr*Str + Gtr;
        V0  = wtp.volume * waterMolarMass;
        Cp0 = wtp.cp * waterMolarMass;
        Cv0 = wtp.cv * waterMolarMass;
        // S0  = Sw + Str;
        // U0  = Uw + Utr;
        // A0  = Uw - T * (Sw + Str) + Ttr * Str + Atr;
    };

    return StandardThermoModel(evalfn, extractParams(params));
}

} // namespace Reaktoro
