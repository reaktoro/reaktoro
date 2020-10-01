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

#include "PhreeqcThermo.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Extensions/Phreeqc/PhreeqcUtils.hpp>
#include <Reaktoro/Thermodynamics/Reactions/ReactionThermoModelAnalyticalPHREEQC.hpp>
#include <Reaktoro/Thermodynamics/Reactions/ReactionThermoModelVantHoff.hpp>

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

    if(species->name == "H2O")
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

auto standardVolume(const PhreeqcPhase* phase, real T, real P) -> real
{
    const auto R = universalGasConstant;
    const auto m3_to_cm3 = 1.0e6;

    if(isGaseousSpecies(phase))
        return R*T/P * m3_to_cm3; // ideal gas volume in cm3/mol

    return phase->logk[vm0]; // constant solid volume in cm3/mol
}

auto pressureCorrectionEnergy(const PhreeqcSpecies* species, real P, real V0) -> real
{
    const auto P0 = 101325.0; // reference pressure in Pa (1 atm = 101325 Pa)
    const auto deltaP = P - P0;
    return (deltaP > 0.0) ? deltaP*V0 : real{0.0}; // correction energy in J/mol
}

auto pressureCorrectionEnergy(const PhreeqcPhase* phase, real P, real V0) -> real
{
    if(isGaseousSpecies(phase))
        return 0.0; // zero energy correction for gases as done in PHREEQC

    const auto P0 = 101325.0; // reference pressure in Pa (1 atm = 101325 Pa)
    const auto deltaP = P - P0;
    return (deltaP > 0.0) ? deltaP*V0 : real{0.0}; // correction energy in J/mol
}

/// Create the standard thermodynamic model of the formation reaction.
template<typename SpeciesType>
auto reactionThermoPropsFnAux(const SpeciesType* s) -> ReactionThermoPropsFn
{
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
        const auto A1 = logk[T_A1];
        const auto A2 = logk[T_A2];
        const auto A3 = logk[T_A3];
        const auto A4 = logk[T_A4];
        const auto A5 = logk[T_A5];
        const auto A6 = logk[T_A6];
        return ReactionThermoModelAnalyticalPHREEQC(A1, A2, A3, A4, A5, A6);
    }
    else
    {
        const auto lgK0 = logk[logK_T0];
        const auto dH0 = logk[delta_h] * 1e3; // convert from kJ/mol to J/mol
        const auto Tref = 298.15; // reference temperature (in K)
        return ReactionThermoModelVantHoff(lgK0, dH0, Tref);
    }
}

auto reactionThermoPropsFn(const PhreeqcSpecies* species) -> ReactionThermoPropsFn
{
    return reactionThermoPropsFnAux(species);
}

auto reactionThermoPropsFn(const PhreeqcPhase* phase) -> ReactionThermoPropsFn
{
    // Note: PHREEQC is not consisent with the direction of the reactions. For
    // gases and minerals, we need to invert the sign of the delta properties
    // of the reaction.
    const auto fn = reactionThermoPropsFnAux(phase);
    return [=](real T, real P)
    {
        ReactionThermoProps props = fn(T, P);
        props.dG0 *= -1.0;
        props.dH0 *= -1.0;
        return props;
    };
}

} // namespace PhreeqcUtils
} // namespace Reaktoro
