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

#include "StandardThermoModelExtendedUNIQUAC.hpp"

// C++ includes
#include <cmath>

// Reaktoro includes
#include <Reaktoro/Serialization/Models/StandardThermoModels.hpp>

namespace Reaktoro {

using std::exp;
using std::log;
using std::pow;
using std::sqrt;

/// Return a Vec<Param> object containing all Param objects in @p params.
auto extractParams(StandardThermoModelParamsExtendedUNIQUAC const& params) -> Vec<Param>
{
    auto const& [Gr, Hr, Sr, Vr, Cp, a, b, c, alpha, beta, Theta] = params;
    return {Gr, Hr, Sr, Vr, Cp, a, b, c, alpha, beta, Theta};
}

/// Return a ModelSerializer for given model parameters in @p params.
auto createModelSerializer(StandardThermoModelParamsExtendedUNIQUAC const& params) -> ModelSerializer
{
    return [=]()
    {
        Data node;
        node["ExtendedUNIQUAC"] = params;
        return node;
    };
}

auto StandardThermoModelExtendedUNIQUAC(StandardThermoModelParamsExtendedUNIQUAC const& params) -> StandardThermoModel
{
    auto evalfn = [=](StandardThermoProps& props, real T, real P)
    {
        // Unpack the model parameters
        auto const& [Gr, Hr, Sr, Vr, Cp, a, b, c, alpha, beta, Theta] = params;

        // Auxiliary variables
        const auto PA_TO_BAR = 1.0e-5;
        const auto KJ_TO_J = 1.0e3;
        const auto Tr    = 298.15;   // referente temperature in K (equivalent to 25 °C)
        const auto Pr    = 1.0;      // referente pressure in bar
        const auto Pb    = P * PA_TO_BAR; // converted pressure from Pa to bar
        const auto T2    = T*T;
        const auto Tr2   = Tr*Tr;
        const auto TmTr  = T - Tr;
        const auto PmPr  = Pb - Pr;
        const auto TmTr2 = TmTr*TmTr;
        const auto PmPr2 = PmPr*PmPr;

        // Compute the standard heat capacity (constant pressure) using the model Cp°(T) = a + bT + c/(T - Θ)
        const auto Cp0 = a + b*T + c/(T - Theta); // in J/(mol·K)

        // Compute the standard enthalpy of the substance using the previous heat capacity model
        const auto H0 = Hr + a*TmTr + 0.5*b*(T2 - Tr2) + c*log((T - Theta)/(Tr - Theta)); // in kJ/(mol·K)

        // Compute the standard Gibbs energy of the substance using the previous heat capacity model
        auto G0 = Gr*T/Tr + Hr*(1 - T/Tr) - a*T*(log(T/Tr) + Tr/T - 1) - 0.5*b*TmTr2 - c*T/Theta*((T - Theta)/T*log((T - Theta)/(Tr - Theta)) - log(T/Tr)); // in kJ/(mol·K)

        // Compute the pressure-correction for the standard Gibbs energy: α(P - Pr) + β(P - Pr)²
        G0 += (alpha + beta*PmPr)*PmPr; // in kJ/(mol·K)

        // Set the corresponding properties in the StandardThermoProps object
        props.Cp0 = Cp0;
        props.H0  = H0 * KJ_TO_J;
        props.G0  = G0 * KJ_TO_J;
    };

    return StandardThermoModel(evalfn, extractParams(params), createModelSerializer(params));
}

} // namespace Reaktoro
