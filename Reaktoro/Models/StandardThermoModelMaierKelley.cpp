// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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

#include "StandardThermoModelMaierKelley.hpp"

// C++ includes
#include <cmath>
using std::log;

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Serialization/Models.YAML.hpp>

namespace Reaktoro {

/// Return a Vec<Param> object containing all Param objects in @p params.
auto extractParams(const StandardThermoModelParamsMaierKelley& params) -> Vec<Param>
{
    const auto& [Gf, Hf, Sr, Vr, a, b, c, Tmax] = params;
    return {Gf, Hf, Sr, Vr, a, b, c};
}

/// Return a ModelSerializer for given model parameters in @p params.
auto createModelSerializer(const StandardThermoModelParamsMaierKelley& params) -> ModelSerializer
{
    return [=]()
    {
        yaml node;
        node["MaierKelley"] = params;
        return node;
    };
}

auto StandardThermoModelMaierKelley(const StandardThermoModelParamsMaierKelley& params) -> StandardThermoModel
{
    const auto isgas = params.Vr.value() == 0.0;

    auto evalfn = [=](StandardThermoProps& props, real T, real P)
    {
        auto& [G0, H0, V0, Cp0, VT0, VP0] = props;
        const auto& [Gf, Hf, Sr, Vr, a, b, c, Tmax] = params;

        const auto Tr = 298.15; // the reference temperature of 25 C (in K)
        const auto Pr = 1.0e5;  // the reference pressure of 1 bar (in Pa)
        const auto R  = universalGasConstant;

        const auto CpdT   = a*(T - Tr) + 0.5*b*(T*T - Tr*Tr) - c*(1.0/T - 1.0/Tr);
        const auto CpdlnT = a*log(T/Tr) + b*(T - Tr) - 0.5*c*(1.0/(T*T) - 1.0/(Tr*Tr));
        const auto VdP    = Vr*(P - Pr);

        V0  = Vr;
        G0  = Gf - Sr*(T - Tr) + CpdT - T*CpdlnT + VdP;
        H0  = Hf + CpdT + VdP;
        Cp0 = a + b*T + c/(T*T);
        VT0 = 0.0;
        VP0 = 0.0;
        // S0  = Sr + CpdlnT;
    };

    return StandardThermoModel(evalfn, extractParams(params), createModelSerializer(params));
}

} // namespace Reaktoro
