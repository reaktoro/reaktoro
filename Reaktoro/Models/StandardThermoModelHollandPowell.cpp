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

#include "StandardThermoModelHollandPowell.hpp"

// C++ includes
#include <cmath>

// Reaktoro includes
#include <Reaktoro/Serialization/SerializationYAML.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterThermoState.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterThermoStateUtils.hpp>

namespace Reaktoro {

using std::exp;
using std::log;
using std::pow;
using std::sqrt;

/// Return a Params object containing all Param objects in @p params.
auto extractParams(const StandardThermoModelParamsHollandPowell& params) -> Params
{
    const auto& [Gf, Hf, Sr, Vr, a, b, c, d, alpha0, kappa0, kappa0p, kappa0pp, numatoms, Tcr, Smax, Vmax, Tmax] = params;
    return {Gf, Hf, Sr, Vr, a, b, c, d, alpha0, kappa0, kappa0p, kappa0pp};
}

/// Return a ModelSerializer for given model parameters in @p params.
auto createModelSerializer(const StandardThermoModelParamsHollandPowell& params) -> ModelSerializer
{
    return [=]()
    {
        yaml node;
        node["HollandPowell"] = params;
        return node;
    };
}

auto StandardThermoModelHollandPowell(const StandardThermoModelParamsHollandPowell& params) -> StandardThermoModel
{
    auto evalfn = [=](StandardThermoProps& props, real T, real P)
    {
        // Unpack the properties to be computed by this model
        auto& [G0, H0, V0, Cp0, Cv0] = props;

        // Unpack the model parameters
        const auto& [Gf, Hf, Sr, Vr, MKa, MKb, MKc, MKd, alpha0, kappa0, kappa0p, kappa0pp, numatoms, Tcr, Smax, Vmax, Tmax] = params;

        // Auxiliary variables related to reference pressure Pr and reference temperature Tr
        const auto Pr   = 1.0e5;
        const auto Tr   = 298.15;
        const auto Tr2  = Tr*Tr;
        const auto Tr05 = sqrt(Tr);

        // Auxiliary variables related to given temperature T
        const auto T2  = T*T;
        const auto T05 = sqrt(T);

        // Compute the heat capacity and the integrals Cp*dT and Cp*d(lnT) according to equations in Fig. 4 of SUPCRTBL (2016)
        const auto Cp     = MKa + MKb*T + MKc/T2 + MKd/T05;
        const auto CpdT   = MKa*(T - Tr) + 0.5*MKb*(T2 - Tr2) - MKc*(1.0/T - 1.0/Tr) + 2.0*MKd*(T05 - Tr05);
        const auto CpdlnT = MKa*log(T/Tr) + MKb*(T - Tr) - 0.5*MKc*(1/T2 - 1/Tr2) - 2.0*MKd*(1/T05 - 1/Tr05);

        // The volume of the substance to be calculated below
        real V = 0.0;

        // The integral(V*dP) to be calculated below
        real VdP = 0.0;

        // Check special case when kappa0 is zero (simpler case)
        if(kappa0 == 0.0)
        {
            // Assume constant volume for the substance
            V = Vr;

            // Compute the integral(V*dP) assuming constant volume
            VdP = Vr*(P - Pr);
        }
        else
        {
            // See equation (3) in Holland and Powell (2011) for the ka, kb, kc (different name to avoid confusion with Maier-Kelley coefficients)
            const auto ka = (1 + kappa0p)/(1 + kappa0p + kappa0*kappa0pp);
            const auto kb = kappa0p/kappa0 - kappa0pp/(1 + kappa0p);
            const auto kc = (1 + kappa0p + kappa0*kappa0pp)/(kappa0p*(1 + kappa0p) - kappa0*kappa0pp);

            // See Holland and Powell (2011) p. 346, end of left column for θ. Sr must be in J/(mol*K)!
            const auto theta  = 10636.0/(Sr/numatoms + 6.44);

            // Define u = θ/T and u0 = θ/Tr (see equation below Fig 1 in Holland and Powell 2011)
            const auto u  = theta/T;
            const auto u0 = theta/Tr;

            // Compute exp(u) and exp(u0) as they are used often
            const auto exp_u  = exp(u);
            const auto exp_u0 = exp(u0);

            // Define w = u/(exp(u) - 1) to simplify computation of ξ (see equation below Fig 1 in Holland and Powell 2011)
            const auto w  = u/(exp_u - 1);
            const auto w0 = u0/(exp_u0 - 1);

            // Compute ξ and ξ0 (see equation below Fig 1 in Holland and Powell 2011)
            const auto E  = w*w*exp_u;
            const auto E0 = w0*w0*exp_u0;

            // Compute Pth value according to the equation (without number) before equation (12) in Holland and Powell (2011)
            const auto Pth = alpha0*kappa0*(theta/E0)*(1/(exp_u - 1) - 1/(exp_u0 - 1));

            // Compute the auxiliary variables below for a more efficient model evaluation
            const auto aux1 = 1 - kb*Pth;
            const auto aux2 = pow(aux1, 1 - kc); // === pow(1 - kb*Pth, 1 - kc)
            const auto aux3 = 1 + kb*(P - Pth);
            const auto aux4 = pow(aux3, 1 - kc); // === pow(1 + kb*(P - Pth), 1 - kc)
            const auto aux5 = kb*(kc - 1)*P;
            const auto aux6 = pow(aux1, kc); // === pow(1 - kb*Pth, kc)        TODO: should be equivalent to aux1/aux2 (avoiding thus an extra exp)
            const auto aux7 = pow(aux3, kc); // === pow(1 + kb*(P - Pth), kc)  TODO: should be equivalent to aux3/aux4 (avoiding thus an extra exp)

            // Compute the thermal expansion α using equation below (12)
            const auto alpha = alpha0*(E/E0)/(aux1*(ka + (1 - ka)*aux6));

            // Compute the compressibility κ using equation below (12)
            const auto kappa = kappa0*aux3*(ka + (1 - ka)*aux7);

            // Compute the volume of the substance according equation (12) in Holland and Powell (2011)
            V = Vr*(1 - ka*(1 - 1/aux7));

            // Compute the integral(V*dP) according to equation (13) in Holland and Powell (2011)
            VdP = P*Vr * (1 - ka + ka*(aux2 - aux4)/aux5);
        }

        // Compute the standard properties unpacked from props
        G0 = Gf - Sr*(T - Tr) + CpdT - T*CpdlnT + VdP; // see equation (2) of SUPCRTBL (2016)
        H0 = Hf + CpdT + VdP; // similar to Maier-Kelley
        V0 = V;
        Cp0 = Cp;
        Cv0 = Cp;
        // S0 = STrPr + CpdlnT;
    };

    return StandardThermoModel(evalfn, extractParams(params), createModelSerializer(params));
}

} // namespace Reaktoro
