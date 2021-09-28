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

#include "StandardThermoModelMineralHKF.hpp"

// C++ includes
#include <cmath>
using std::log;

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Serialization/Models.YAML.hpp>

namespace Reaktoro {

/// Return a Vec<Param> object containing all Param objects in @p params.
auto extractParams(const StandardThermoModelParamsMineralHKF& params) -> Vec<Param>
{
    const auto& [Gf, Hf, Sr, Vr, ntr, a, b, c, Ttr, Htr, Vtr, dPdTtr, Tmax] = params;
    Vec<Param> collected = {Gf, Hf, Sr, Vr};
    for(const auto& x : params.a) collected.push_back(x);
    for(const auto& x : params.b) collected.push_back(x);
    for(const auto& x : params.c) collected.push_back(x);
    for(const auto& x : params.Ttr) collected.push_back(x);
    for(const auto& x : params.Htr) collected.push_back(x);
    for(const auto& x : params.Vtr) collected.push_back(x);
    for(const auto& x : params.dPdTtr) collected.push_back(x);
    return collected; // TODO: Implement constructor Vec<Param>(const Args&... args) so that we can write `return Vec<Param>(Gf, Hf, Sr, Vr, a, b, c, Ttr, Htr, Vtr, dPdTtr);`
}

/// Return a ModelSerializer for given model parameters in @p params.
auto createModelSerializer(const StandardThermoModelParamsMineralHKF& params) -> ModelSerializer
{
    return [=]()
    {
        yaml node;
        node["MineralHKF"] = params;
        return node;
    };
}

auto StandardThermoModelMineralHKF(const StandardThermoModelParamsMineralHKF& params) -> StandardThermoModel
{
    auto evalfn = [=](StandardThermoProps& props, real T, real P)
    {
        auto& [G0, H0, V0, Cp0, VT0, VP0] = props;
        const auto& [Gf, Hf, Sr, Vr, ntr, a, b, c, Ttr, Htr, Vtr, dPdTtr, Tmax] = params;

        // Auxiliary variables
        const auto Tr = 298.15;  // The reference temperature (in K)
        const auto Pr = 1.0e+05; // The reference pressure (in Pa)

        // Collect the temperature points used for the integrals along the pressure line P = Pr
        Vec<real> Ti;
        Ti.push_back(Tr);
        for(int i = 0; i < ntr; ++i)
            if(T > Ttr[i]) Ti.push_back(Ttr[i]);
        Ti.push_back(T);

        // Collect the pressure intercepts along the temperature line T for every phase transition boundary
        Vec<real> Ptr;
        for(int i = 0; i < ntr; ++i)
        {
            if(dPdTtr[i] != 0.0)
                Ptr.push_back(Pr + dPdTtr[i]*(T - Ttr[i]));
        }

        // Calculate the heat capacity of the mineral at T
        real Cp = 0.0;
        for(auto i = 0; i+1 < Ti.size(); ++i)
            if(Ti[i] <= T && T <= Ti[i+1])
                Cp = a[i] + b[i]*T + c[i]/(T*T); // Cp evaluated using coefficients of the appropriate interval (see `SUBROUTINE Cptrms` in supcrt92/reac92d.f)

        // Calculate the integrals of the heat capacity function of the mineral from Tr to T at constant pressure Pr
        real CpdT = 0.0;
        real CpdlnT = 0.0;
        for(auto i = 0; i+1 < Ti.size(); ++i)
        {
            const auto T0 = Ti[i];
            const auto T1 = Ti[i+1];

            CpdT += a[i]*(T1 - T0) + 0.5*b[i]*(T1*T1 - T0*T0) - c[i]*(1.0/T1 - 1.0/T0); // see `FUNCTION CpdT` in supcrt92/reac92d.f
            CpdlnT += a[i]*log(T1/T0) + b[i]*(T1 - T0) - 0.5*c[i]*(1.0/(T1*T1) - 1.0/(T0*T0)); // see `FUNCTION CpdlnT` in supcrt92/reac92d.f
        }

        // Calculate the volume and other auxiliary quantities for the thermodynamic properties of the mineral
        real V = Vr;
        real GdH = 0.0; // last term in equation (82) of SUPCRT92 paper
        real HdH = 0.0; // last term in equation (79) of SUPCRT92 paper
        real SdH = 0.0; // last term in equation (80) of SUPCRT92 paper
        for(unsigned i = 1; i+1 < Ti.size(); ++i)
        {
            GdH += Htr[i-1]*(T - Ti[i])/Ti[i]; // see `SUBROUTINE pttrms` in supcrt92/reac92d.f
            HdH += Htr[i-1];
            SdH += Htr[i-1]/Ti[i];

            V += Vtr[i-1];
        }

        // Calculate the volume integral from Pr to P at constant temperature T
        real VdP = V*(P - Pr); // start with full VdP and decrease accordingly below due to phase transitions
        for(unsigned i = 0; i < Ptr.size(); ++i)
        {
            if(0.0 < Ptr[i] && Ptr[i] < P)
            {
                V   -= Vtr[i];
                VdP -= Vtr[i]*(P - Ptr[i]);
            }
        }

        // Calculate the standard molal thermodynamic properties of the mineral
        V0 = V;
        G0 = Gf - Sr*(T - Tr) + CpdT - T*CpdlnT + VdP - GdH;
        H0 = Hf + CpdT + VdP + HdH;
        Cp0 = Cp;
        VT0 = 0.0;
        VP0 = 0.0;
        // S0 = Sr + CpdlnT + SdH;
    };

    return StandardThermoModel(evalfn, extractParams(params), createModelSerializer(params));
}

} // namespace Reaktoro
