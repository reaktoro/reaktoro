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

#include "PhreeqcThermo.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Extensions/Phreeqc/PhreeqcUtils.hpp>
#include <Reaktoro/Extensions/Phreeqc/PhreeqcWater.hpp>
#include <Reaktoro/Models/StandardThermoModels/ReactionStandardThermoModelConstLgK.hpp>
#include <Reaktoro/Models/StandardThermoModels/ReactionStandardThermoModelPhreeqcLgK.hpp>
#include <Reaktoro/Models/StandardThermoModels/ReactionStandardThermoModelVantHoff.hpp>

namespace Reaktoro {
namespace PhreeqcUtils {

auto standardVolume(const PhreeqcSpecies* species, real T, real P, const PhreeqcWaterProps& wprops) -> real
{
    //--------------------------------------------------------------------------------
    // Implementation based on PHREEQC method `Phreeqc::calc_vm`
    //--------------------------------------------------------------------------------

    // COMMENT FROM PHREEQC
    //
    // Calculate molar volumes for aqueous species with a Redlich type eqn:
    //  Vm = Vm0(tc) + (Av / 2) * z^2 * I^0.5 + coef(tc) * I^(b4).
    //   Vm0(tc) is calc'd using supcrt parms, or from millero[0] + millero[1] * tc + millero[2] * tc^2
    //   for Av * z^2 * I^0.5, see Redlich and Meyer, Chem. Rev. 64, 221.
    //      Av is in (cm3/mol)(mol/kg)^-0.5, = DH_Av.
    //      If b_Av != 0, the extended DH formula is used: I^0.5 /(1 + b_Av * DH_B * I^0.5).
    //      DH_Av and DH_B are from calc_dielectrics(tc, pa).
    // coef(tc) = logk[vmi1] + logk[vmi2] / (TK - 228) + logk[vmi3] * (TK - 228).
    //   b4 = logk[vmi4], or
    // coef(tc) = millero[3] + millero[4] * tc + millero[5] * tc^2

    const auto& [wtp, wep] = wprops;
    const auto& [rho_0, kappa_0] = wtp; // references to the given thermo properties of water
    const auto& [eps_r, DH_A, DH_B, DH_Av, ZBrn, QBrn] = wep; // references to given electro properties of water

    const auto pascal_to_atm = 9.86923e-6;
    const auto cm3_to_m3 = 1e-6;

    const auto tc = T - 273.15;
    const auto pa = P * pascal_to_atm;

    const auto pb_s = 2600. + pa * 1.01325;
    const auto TK_s = tc + 45.15;

    // Note: In Reaktoro, standard thermo properties do not depend on
    // concentration variables. There is thus no dependence on ionic strength
    // in the calculation of standard molar volumes of the species below.

    const auto a1 = species->logk[vma1];
    const auto a2 = species->logk[vma2];
    const auto a3 = species->logk[vma3];
    const auto a4 = species->logk[vma4];
    const auto wr = species->logk[wref];

    if(PhreeqcUtils::name(species) == "H2O")
        return 18.016 / rho_0; // in cm3/mol

    if(species->logk[vma1])
        return a1 + a2 / pb_s + (a3 + a4 / pb_s) / TK_s - wr*QBrn; // in cm3/mol

    if(species->millero[0])
    {
        const auto m0 = species->millero[0];
        const auto m1 = species->millero[1];
        const auto m2 = species->millero[2];
        return m0 + tc*(m1 + tc*m2); // in cm3/mol
    }

    return 0.0;
}

auto standardVolumeIonicStrengthCorrection(const PhreeqcSpecies* species, real T, real P, real mu, const PhreeqcWaterProps& wprops) -> real
{
    //--------------------------------------------------------------------------------
    // Implementation based on PHREEQC method `Phreeqc::calc_vm`
    //--------------------------------------------------------------------------------

    const auto z = PhreeqcUtils::charge(species);

    if(z == 0.0)
        return 0.0;

    const auto sqrt_mu = sqrt(mu);

    const auto& DH_B = wprops.wep.DH_B;
    const auto& DH_Av = wprops.wep.DH_Av;

    const auto tc = T - 273.15;
    const auto TK_s = tc + 45.15;

    const auto i1 = species->logk[vmi1];
    const auto i2 = species->logk[vmi2];
    const auto i3 = species->logk[vmi3];
    const auto i4 = species->logk[vmi4];

    if(species->logk[vma1])
    {
        real corr = 0.0;

        const auto bAv = species->logk[b_Av];

        if(bAv < 1e-5)
            corr += z * z * 0.5 * DH_Av * sqrt_mu;
        else
            corr += z * z * 0.5 * DH_Av * sqrt_mu / (1 + bAv * DH_B * sqrt_mu);

        if(i1 != 0.0 || i2 != 0.0 || i3 != 0.0)
        {
            real bi = i1 + i2 / TK_s + i3 * TK_s;
            corr += i4 == 1.0 ? bi * mu : bi * pow(mu, i4);
        }

        return corr;
    }

    if(species->millero[0])
    {
        const auto m3 = species->millero[3];
        const auto m4 = species->millero[4];
        const auto m5 = species->millero[5];

        return z * z * 0.5 * DH_Av * sqrt_mu + (m3 + tc * (m4 + tc * m5)) * mu;
    }

    return 0.0;
}

auto standardVolume(const PhreeqcPhase* phase, real T, real P) -> real
{
    return phase->logk[vm0]; // constant solid volume in cm3/mol or zero volume for gases
}

/// Create the standard thermodynamic model of the formation reaction.
template<typename SpeciesType>
auto reactionThermoModelAux(const SpeciesType* s, double sign) -> ReactionStandardThermoModel
{
    if(PhreeqcUtils::reactants(s).empty())
    {
        ReactionStandardThermoModelParamsConstLgK params;
        params.lgKr = 0.0;
        params.Pr = 101'325; // reference pressure (1 atm = 101325 Pa)
        return ReactionStandardThermoModelConstLgK(params);
    }

    const auto logk = s->logk;

    const auto use_analytic_expression = [&]() -> bool
    {
        for(int i = T_A1; i <= T_A6; ++i)
            if(logk[i] != 0.0)
                return true;
        return false;
    };

    if(use_analytic_expression())
    {
        ReactionStandardThermoModelParamsPhreeqcLgK params;
        params.A1 = sign * logk[T_A1];
        params.A2 = sign * logk[T_A2];
        params.A3 = sign * logk[T_A3];
        params.A4 = sign * logk[T_A4];
        params.A5 = sign * logk[T_A5];
        params.A6 = sign * logk[T_A6];
        params.Pr = 101'325; // reference pressure (1 atm = 101325 Pa)
        return ReactionStandardThermoModelPhreeqcLgK(params);
    }
    else
    {
        ReactionStandardThermoModelParamsVantHoff params;
        params.lgKr = sign * logk[logK_T0];
        params.dHr = sign * logk[delta_h] * 1e3; // convert from kJ/mol to J/mol
        params.Tr = 298.15; // reference temperature (in K)
        params.Pr = 101'325; // reference pressure (1 atm = 101325 Pa)
        return ReactionStandardThermoModelVantHoff(params);
    }
}

auto reactionThermoModel(const PhreeqcSpecies* species) -> ReactionStandardThermoModel
{
    const auto sign = 1.0;
    return reactionThermoModelAux(species, sign);
}

auto reactionThermoModel(const PhreeqcPhase* phase) -> ReactionStandardThermoModel
{
    // Note: PHREEQC is not consisent with the direction of the reactions. For
    // gases and minerals, we need to invert the sign of the delta properties
    // of the reaction.
    const auto sign = -1.0;
    return reactionThermoModelAux(phase, sign);
}

auto standardVolumeModel(const PhreeqcSpecies* species) -> Model<real(real,real)>
{
    return Model<real(real,real)>([=](real T, real P) -> real
    {
        const auto wprops = waterPropsMemoized(T, P);
        return standardVolume(species, T, P, wprops) * cubicCentimeterToCubicMeter;
    });
}

auto standardVolumeModel(const PhreeqcPhase* phase) -> Model<real(real,real)>
{
    return Model<real(real,real)>([=](real T, real P) -> real
    {
        return standardVolume(phase, T, P) * cubicCentimeterToCubicMeter;
    });
}

} // namespace PhreeqcUtils
} // namespace Reaktoro
