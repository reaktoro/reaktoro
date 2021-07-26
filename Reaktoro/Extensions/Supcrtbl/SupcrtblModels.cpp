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

#include "SupcrtblModels.hpp"

// C++ includes
#include <cmath>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/ConvertUtils.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/NamingUtils.hpp>
#include <Reaktoro/Extensions/Supcrt/SpeciesElectroProps.hpp>
#include <Reaktoro/Extensions/Supcrt/SpeciesElectroPropsHKF.hpp>
#include <Reaktoro/Extensions/Supcrt/SpeciesThermoState.hpp>
#include <Reaktoro/Extensions/Supcrt/SupcrtModels.hpp>
#include <Reaktoro/Extensions/Supcrtbl/SupcrtblParams.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterElectroProps.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterElectroPropsJohnsonNorton.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterThermoState.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterThermoStateUtils.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterUtils.hpp>

namespace Reaktoro {

using std::log;
using std::sqrt;

namespace {

const auto cal_to_J  = 4.184;
const auto cal_to_kJ = 4.184e-3;
const auto kJ_to_cal = 1/cal_to_kJ;
const auto kJ_to_J   = 1e+3;
const auto J_to_kJ   = 1e-3;
const auto Pa_to_bar = 1e-5;

const auto Tr   = 298.15;
const auto Tr2  = Tr*Tr;
const auto Tr05 = sqrt(Tr);

} // namespace


auto supcrtblStandardThermoPropsSolventHKF(real T, real P, const SupcrtblParamsAqueousSolventHKF& params, const WaterThermoState& wts) -> SpeciesThermoState
{
    SupcrtParamsAqueousSolventHKF cparams;
    // cparams.Ttr = params.Ttr; // TODO: Revisit this because currently SupcrtParamsAqueousSolventHKF is empty (no data members)
    // cparams.Str = params.Str;
    // cparams.Gtr = params.Gtr;
    // cparams.Htr = params.Htr;
    // cparams.Utr = params.Utr;
    // cparams.Atr = params.Atr;
    return supcrtStandardThermoPropsSolventHKF(T, P, cparams, wts);
}

auto supcrtblStandardThermoPropsSoluteHKF(real T, real P, const SupcrtblParamsAqueousSoluteHKF& params, const SpeciesElectroProps& aes, const WaterElectroProps& wes) -> SpeciesThermoState
{
    SupcrtParamsAqueousSoluteHKF cparams;
    cparams.name = params.name;
    cparams.charge = params.charge;
    cparams.Gf = params.Gf * kJ_to_cal;
    cparams.Hf = params.Hf * kJ_to_cal;
    cparams.Sr = params.Sr;
    cparams.a1 = params.a1 * kJ_to_cal * 1e-1;
    cparams.a2 = params.a2 * kJ_to_cal * 1e+2;
    cparams.a3 = params.a3 * kJ_to_cal;
    cparams.a4 = params.a4 * kJ_to_cal * 1e+4;
    cparams.c1 = params.c1 * kJ_to_cal;
    cparams.c2 = params.c2 * kJ_to_cal * 1e+4;
    cparams.wref = params.wref * 1e+5;
    return supcrtStandardThermoPropsSoluteHKF(T, P, cparams, aes, wes);
}

auto supcrtblStandardThermoPropsFluidHollandPowell(real T, real P, const SupcrtblParamsFluidHollandPowell& params) -> SpeciesThermoState
{
    const auto& [Gf, Hf, Sr, Vr, a, b, c, d, Tmax] = params;

    const auto STrPr = Sr * J_to_kJ;

    const auto T2  = T*T;
    const auto T05 = sqrt(T);

    const auto Cp     = a + b*T + c/T2 + d/T05;
    const auto CpdT   = a*(T - Tr) + b/2.0*(T2 - Tr2) - c*(1.0/T - 1.0/Tr) + d/0.5*(T05 - Tr05);
    const auto CpdlnT = a*log(T/Tr) + b*(T - Tr) - c/2.0*(1/T2 - 1/Tr2) - d/0.5*(1/T05 - 1/Tr05);

    const auto R = universalGasConstant;

    auto S  = STrPr + CpdlnT;
    auto H  = Hf + CpdT;
    auto G  = Gf - STrPr*(T - Tr) + CpdT - T*CpdlnT;
    auto U  = H;
    auto A  = U - T*S;

    S *= kJ_to_J;
    H *= kJ_to_J;
    G *= kJ_to_J;
    U *= kJ_to_J;
    A *= kJ_to_J;

    SpeciesThermoState props;
    props.volume           = 0.0;
    props.gibbs_energy     = G;
    props.enthalpy         = H;
    props.entropy          = S;
    props.internal_energy  = U;
    props.helmholtz_energy = A;
    props.heat_capacity_cp = Cp;
    props.heat_capacity_cv = props.heat_capacity_cp - R;

    return props;
}

auto supcrtblStandardThermoPropsMineralHollandPowell(real T, real P, const SupcrtblParamsMineralHollandPowell& params) -> SpeciesThermoState
{
    const auto& [Gf, Hf, Sr, Vr, MKa, MKb, MKc, MKd, alpha0, kappa0, kappap0, kappapp0, numatoms, Tmax] = params;

    const auto STrPr = Sr * J_to_kJ;

    const auto T2  = T*T;
    const auto T05 = sqrt(T);
    const auto Pbar = P * Pa_to_bar;

    const auto ka     = (1 + kappap0)/(1 + kappap0 + kappa0*kappapp0);
    const auto kb     = kappap0/kappa0 - kappapp0/(1 + kappap0);
    const auto kc     = (1 + kappap0 + kappa0*kappapp0)/(kappap0*(1 + kappap0) - kappa0*kappapp0);
    const auto theta  = 10636.0/(Sr/numatoms + 6.44); // Sr must be in J/(mol*K)!
    const auto u      = theta/T;
    const auto u0     = theta/Tr;
    const auto exp_u  = exp(u);
    const auto exp_u0 = exp(u0);
    const auto w      = u/(exp_u - 1);
    const auto E      = w*w*exp_u;
    const auto w0     = u0/(exp_u0 - 1);
    const auto E0     = w0*w0*exp_u0;
    const auto Pth    = alpha0*kappa0*(theta/E0)*(1/(exp_u - 1) - 1/(exp_u0 - 1)); // TODO: Note kappa0 is in kbar (that's why SUPCRTBL has a 0.001 factor on P, because their P is in bar)

    const auto aux1 = 1 - kb*Pth;
    const auto aux2 = pow(aux1, 1 - kc); // === pow(1 - kb*Pth, 1 - kc)
    const auto aux3 = 1 + kb*(Pbar - Pth);
    const auto aux4 = pow(aux3, 1 - kc); // === pow(1 + kb*(Pbar - Pth), 1 - kc)
    const auto aux5 = kb*(kc - 1)*Pbar;
    const auto aux6 = pow(aux1, kc); // === pow(1 - kb*Pth, kc)                TODO: should be equivalent to aux1/aux2 (avoiding thus an extra exp)
    const auto aux7 = pow(aux3, kc); // === pow(1 + kb*(Pbar - Pth), kc)       TODO: should be equivalent to aux3/aux4 (avoiding thus an extra exp)

    const auto VdP  = Pbar*Vr * (1 - ka + ka*(aux2 - aux4)/aux5);

    const auto V = Vr*(1 - ka*(1 - 1/aux7));

    const auto alpha = alpha0*(E/E0)/(aux1*(ka + (1 - ka)*aux6));
    const auto kappa = kappa0*aux3*(ka + (1 - ka)*aux7);

    const auto Cp     = MKa + MKb*T + MKc/T2 + MKd/T05;
    const auto CpdT   = MKa*(T - Tr) + MKb/2.0*(T2 - Tr2) - MKc*(1.0/T - 1.0/Tr) + MKd/0.5*(T05 - Tr05);
    const auto CpdlnT = MKa*log(T/Tr) + MKb*(T - Tr) - MKc/2.0*(1/T2 - 1/Tr2) - MKd/0.5*(1/T05 - 1/Tr05);

    const auto S = STrPr + CpdlnT;
    const auto H = Hf + CpdT + VdP;
    const auto G = Gf - STrPr*(T - Tr) + CpdT - T*CpdlnT + VdP;
    const auto U = H - Pbar*V;
    const auto A = U - T*S;

    SpeciesThermoState props;
    props.volume           = V * 1e-5; // from J/(bar*mol) to m3/mol
    props.gibbs_energy     = G * kJ_to_J;
    props.enthalpy         = H * kJ_to_J;
    props.entropy          = S * kJ_to_J;
    props.internal_energy  = U * kJ_to_J;
    props.helmholtz_energy = A * kJ_to_J;
    props.heat_capacity_cp = Cp * kJ_to_J;
    props.heat_capacity_cv = props.heat_capacity_cp; // approximate Cp = Cv for a solid

    return props;
}

// auto supcrtblStandardThermoPropsMineralHollandPowellLandau(real T, real P, const SupcrtblParamsMineralHollandPowellLandau& params) -> SpeciesThermoState
// {

// }


// using std::log;
// using std::min;
// using std::pow;

// namespace {

// /// The reference temperature assumed in the HKF equations of state (in units of K)
// const auto referenceTemperature = 298.15;

// /// The reference temperature assumed in the HKF equations of state (in units of bar)
// const auto referencePressure = 1.0;

// /// The reference dielectric constant of water \epsilon
// const auto referenceDielectricConstant = 78.24385513;

// /// The reference Born function Z (dimensionless)
// const auto referenceBornZ = -1.278055636e-02;

// /// The reference Born function Y (dimensionless)
// const auto referenceBornY = -5.795424563e-05;

// /// The reference Born function Q (dimensionless)
// const auto referenceBornQ =  6.638388994e-12;

// /// The reference Born function N (dimensionless)
// const auto referenceBornN = -2.321814455e-20;

// /// The reference Born function U (dimensionless)
// const auto referenceBornU = 4.872982291e-14;

// /// The reference Born function X (dimensionless)
// const auto referenceBornX = -3.060388224e-07;

// /// The \eta constant in the HKF model (in units of (A*cal)/mol)
// const auto eta = 1.66027e+05;

// /// The constant characteristics \Theta of the solvent (in units of K)
// const auto theta = 228.0;

// /// The constant characteristics \Psi of the solvent (in units of bar)
// const auto psi = 2600.0;

// template<class SupcrtParamsType>
// auto checkTemperatureValidityHKF(real& T, const SupcrtParamsType& params) -> void
// {
//     // Check if temperature bounds should be enforced
//     Assert(T >= 0.0 && T <= params.Tmax, "Unable to calculate the "
//         "thermodynamic properties of species" + params.name + "using the "
//             "revised HKF equations of state.", "The provided temperature `" +
//                 std::to_string(T) + " K` is either negative or greater than the "
//                     "maximum allowed `" + std::to_string(params.Tmax) + " K`.");

//     // Ensure temperature is not above the maximum allowed
//     T = min(T, params.Tmax);
// }

// } // namespace

// auto supcrtStandardThermoPropsSolventHKF(real T, real P, const SupcrtParamsAqueousSolventHKF& params, const WaterThermoState& wt) -> SpeciesThermoState
// {
//     // Auxiliary data from Helgeson and Kirkham (1974), on page 1098
//     const auto Ttr =  273.16;                   // unit: K
//     const auto Str =  15.1320 * calorieToJoule; // unit: J/(mol*K)
//     const auto Gtr = -56290.0 * calorieToJoule; // unit: J/mol
//     const auto Htr = -68767.0 * calorieToJoule; // unit: J/mol
//     const auto Utr = -67887.0 * calorieToJoule; // unit: J/mol
//     const auto Atr = -55415.0 * calorieToJoule; // unit: J/mol

//     const auto Sw = waterMolarMass * wt.entropy;         // unit: J/(mol*K)
//     const auto Hw = waterMolarMass * wt.enthalpy;        // unit: J/mol
//     const auto Uw = waterMolarMass * wt.internal_energy; // unit: J/mol

//     // Calculate the standard molal thermodynamic properties of the aqueous species
//     const auto S  = Sw + Str;
//     const auto H  = Hw + Htr;
//     const auto U  = Uw + Utr;
//     const auto G  = Hw - T * (Sw + Str) + Ttr * Str + Gtr;
//     const auto A  = Uw - T * (Sw + Str) + Ttr * Str + Atr;
//     const auto V  = wt.volume * waterMolarMass;
//     const auto Cp = wt.cp * waterMolarMass;
//     const auto Cv = wt.cv * waterMolarMass;

//     SpeciesThermoState state;
//     state.entropy          = S;
//     state.enthalpy         = H;
//     state.internal_energy  = U;
//     state.gibbs_energy     = G;
//     state.helmholtz_energy = A;
//     state.volume           = V;
//     state.heat_capacity_cp = Cp;
//     state.heat_capacity_cv = Cv;

//     return state;
// }

// auto supcrtStandardThermoPropsSoluteHKF(real T, real P, const SupcrtParamsAqueousSoluteHKF& params, const SpeciesElectroProps& aes, const WaterElectroProps& wes) -> SpeciesThermoState
// {
//     // Auxiliary variables
//     const auto Pbar = P * 1.0e-05;
//     const auto Tr   = referenceTemperature;
//     const auto Pr   = referencePressure;
//     const auto Zr   = referenceBornZ;
//     const auto Yr   = referenceBornY;
//     const auto Gf   = params.Gf;
//     const auto Hf   = params.Hf;
//     const auto Sr   = params.Sr;
//     const auto a1   = params.a1;
//     const auto a2   = params.a2;
//     const auto a3   = params.a3;
//     const auto a4   = params.a4;
//     const auto c1   = params.c1;
//     const auto c2   = params.c2;
//     const auto wr   = params.wref;
//     const auto w    = aes.w;
//     const auto wT   = aes.wT;
//     const auto wP   = aes.wP;
//     const auto wTT  = aes.wTT;
//     const auto Z    = wes.bornZ;
//     const auto Y    = wes.bornY;
//     const auto Q    = wes.bornQ;
//     const auto X    = wes.bornX;

//     // Calculate the standard molal thermodynamic properties of the aqueous species
//     auto V = a1 + a2/(psi + Pbar) +
//         (a3 + a4/(psi + Pbar))/(T - theta) - w*Q - (Z + 1)*wP;

//     auto G = Gf - Sr*(T - Tr) - c1*(T*log(T/Tr) - T + Tr)
//         + a1*(Pbar - Pr) + a2*log((psi + Pbar)/(psi + Pr))
//         - c2*((1.0/(T - theta) - 1.0/(Tr - theta))*(theta - T)/theta
//         - T/(theta*theta)*log(Tr/T * (T - theta)/(Tr - theta)))
//         + 1.0/(T - theta)*(a3*(Pbar - Pr) + a4*log((psi + Pbar)/(psi + Pr)))
//         - w*(Z + 1) + wr*(Zr + 1) + wr*Yr*(T - Tr);

//     auto H = Hf + c1*(T - Tr) - c2*(1.0/(T - theta) - 1.0/(Tr - theta))
//         + a1*(Pbar - Pr) + a2*log((psi + Pbar)/(psi + Pr))
//         + (2.0*T - theta)/pow(T - theta, 2)*(a3*(Pbar - Pr)
//         + a4*log((psi + Pbar)/(psi + Pr)))
//         - w*(Z + 1) + w*T*Y + T*(Z + 1)*wT + wr*(Zr + 1) - wr*Tr*Yr;

//     auto S = Sr + c1*log(T/Tr) - c2/theta*(1.0/(T - theta)
//         - 1.0/(Tr - theta) + log(Tr/T * (T - theta)/(Tr - theta))/theta)
//         + 1.0/pow(T - theta, 2)*(a3*(Pbar - Pr) + a4*log((psi + Pbar)/(psi + Pr)))
//         + w*Y + (Z + 1)*wT - wr*Yr;

//     auto Cp = c1 + c2/pow(T - theta, 2) - (2.0*T/pow(T - theta, 3))*(a3*(Pbar - Pr)
//         + a4*log((psi + Pbar)/(psi + Pr))) + w*T*X + 2.0*T*Y*wT + T*(Z + 1.0)*wTT;

//     auto U = H - Pbar*V;

//     auto A = U - T*S;

//     // Convert the thermodynamic properties of the gas to the standard units
//     V  *= calorieToJoule/barToPascal;
//     G  *= calorieToJoule;
//     H  *= calorieToJoule;
//     S  *= calorieToJoule;
//     U  *= calorieToJoule;
//     A  *= calorieToJoule;
//     Cp *= calorieToJoule;

//     SpeciesThermoState state;
//     state.volume           = V;
//     state.gibbs_energy     = G;
//     state.enthalpy         = H;
//     state.entropy          = S;
//     state.internal_energy  = U;
//     state.helmholtz_energy = A;
//     state.heat_capacity_cp = Cp;
//     state.heat_capacity_cv = state.heat_capacity_cp; // approximate Cp = Cv for an aqueous solution

//     return state;
// }

// // auto supcrtStandardThermoPropsSoluteHKF(real T, real P, const SupcrtParamsAqueousSoluteHKF& params) -> SpeciesThermoState
// // {
// //     WaterThermoState wt = waterThermoStateWagnerPruss(T, P, StateOfMatter::Liquid);

// //     WaterElectroProps wes = waterElectroPropsJohnsonNorton(T, P, wt);

// //     FunctionG g = functionG(T, P, wt);

// //     SpeciesElectroProps aes = speciesElectroPropsHKF(g, params);

// //     return supcrtStandardThermoPropsSoluteHKF(T, P, params, aes, wes);
// // }

// auto supcrtStandardThermoPropsMaierKelly(real T, real P, const SupcrtParamsMaierKelly& params) -> SpeciesThermoState
// {
//     // Check temperature range validity
//     checkTemperatureValidityHKF(T, params);

//     // Auxiliary variables
//     const auto R    = universalGasConstant;
//     const auto Pbar = P * 1.0e-5;
//     const auto Tr   = referenceTemperature;
//     const auto Gf   = params.Gf;
//     const auto Hf   = params.Hf;
//     const auto Sr   = params.Sr;
//     const auto a    = params.a;
//     const auto b    = params.b;
//     const auto c    = params.c;

//     // Calculate the integrals of the heal capacity function of the gas from Tr to T at constant pressure Pr
//     const auto CpdT   = a*(T - Tr) + 0.5*b*(T*T - Tr*Tr) - c*(1.0/T - 1.0/Tr);
//     const auto CpdlnT = a*log(T/Tr) + b*(T - Tr) - 0.5*c*(1.0/(T*T) - 1.0/(Tr*Tr));

//     // Calculate the standard molal thermodynamic properties of the gas
//     auto V  = R*T/P; // the ideal gas molar volume (in units of m3/mol)
//     auto G  = Gf - Sr * (T - Tr) + CpdT - T * CpdlnT;
//     auto H  = Hf + CpdT;
//     auto S  = Sr + CpdlnT;
//     auto U  = H - Pbar*V;
//     auto A  = U - T*S;
//     auto Cp = a + b*T + c/(T*T);

//     // Convert the thermodynamic properties of the gas to the standard units
//     G  *= calorieToJoule;
//     H  *= calorieToJoule;
//     S  *= calorieToJoule;
//     U  *= calorieToJoule;
//     A  *= calorieToJoule;
//     Cp *= calorieToJoule;

//     SpeciesThermoState state;
//     state.volume           = V;
//     state.gibbs_energy     = G;
//     state.enthalpy         = H;
//     state.entropy          = S;
//     state.internal_energy  = U;
//     state.helmholtz_energy = A;
//     state.heat_capacity_cp = Cp;
//     state.heat_capacity_cv = state.heat_capacity_cp - R;

//     return state;
// }

// auto supcrtStandardThermoPropsMaierKellyHKF(real T, real P, const SupcrtParamsMaierKellyHKF& params) -> SpeciesThermoState
// {
//     // Check temperature range validity
//     checkTemperatureValidityHKF(T, params);

//     // Auxiliary variables
//     const auto  Pb   = P * 1.0e-5;
//     const auto& Tr   = referenceTemperature;
//     const auto& Pr   = referencePressure;
//     const auto& Gf   = params.Gf;
//     const auto& Hf   = params.Hf;
//     const auto& Sr   = params.Sr;
//     const auto& Vr   = params.Vr;
//     const auto& nt   = params.nptrans;
//     const auto& a    = params.a;
//     const auto& b    = params.b;
//     const auto& c    = params.c;
//     const auto& Tt   = params.Ttr;
//     const auto& dHt  = params.Htr;
//     const auto& dVt  = params.Vtr;
//     const auto& dPdT = params.dPdTtr;

//     // Collect the temperature points used for the integrals along the pressure line P = Pr
//     std::vector<real> Ti;

//     Ti.push_back(Tr);

//     for(int i = 0; i < nt; ++i)
//         if(T > Tt[i]) Ti.push_back(Tt[i]);

//     Ti.push_back(T);

//     // Collect the pressure intercepts along the temperature line T for every phase transition boundary (see
//     std::vector<real> Pt;
//     for(int i = 0; i < nt; ++i)
//     {
//         if(dPdT[i] != 0.0)
//             Pt.push_back(Pr + dPdT[i]*(T - Tt[i]));
//     }

//     // Calculate the heat capacity of the mineral at T
//     real Cp = {};
//     for(unsigned i = 0; i+1 < Ti.size(); ++i)
//         if(Ti[i] <= T && T <= Ti[i+1])
//             Cp = a[i] + b[i]*T + c[i]/(T*T);

//     // Calculate the integrals of the heat capacity function of the mineral from Tr to T at constant pressure Pr
//     real CpdT = {};
//     real CpdlnT = {};
//     for(unsigned i = 0; i+1 < Ti.size(); ++i)
//     {
//         const auto T0 = Ti[i];
//         const auto T1 = Ti[i+1];

//         CpdT += a[i]*(T1 - T0) + 0.5*b[i]*(T1*T1 - T0*T0) - c[i]*(1.0/T1 - 1.0/T0);
//         CpdlnT += a[i]*log(T1/T0) + b[i]*(T1 - T0) - 0.5*c[i]*(1.0/(T1*T1) - 1.0/(T0*T0));
//     }

//     // Calculate the volume and other auxiliary quantities for the thermodynamic properties of the mineral
//     real V = Vr;
//     real GdH = {};
//     real HdH = {};
//     real SdH = {};
//     for(unsigned i = 1; i+1 < Ti.size(); ++i)
//     {
//         GdH += dHt[i-1]*(T - Ti[i])/Ti[i];
//         HdH += dHt[i-1];
//         SdH += dHt[i-1]/Ti[i];

//         V += dVt[i-1];
//     }

//     // Calculate the volume integral from Pr to P at constant temperature T
//     real VdP = 0.023901488*V*(Pb - Pr);
//     for(unsigned i = 0; i < Pt.size(); ++i)
//     {
//         if(0.0 < Pt[i] && Pt[i] < Pb)
//         {
//             V   -= dVt[i];
//             VdP -= 0.023901488*dVt[i]*(Pb - Pt[i]);
//         }
//     }

//     // Calculate the standard molal thermodynamic properties of the mineral
//     auto G = Gf - Sr * (T - Tr) + CpdT - T * CpdlnT + VdP - GdH;
//     auto H = Hf + CpdT + VdP + HdH;
//     auto S = Sr + CpdlnT + SdH;
//     auto U = H - Pb*V;
//     auto A = U - T*S;

//     // Convert the thermodynamic properties of the mineral to the standard
//     V  *= cubicCentimeterToCubicMeter;
//     G  *= calorieToJoule;
//     H  *= calorieToJoule;
//     S  *= calorieToJoule;
//     U  *= calorieToJoule;
//     A  *= calorieToJoule;
//     Cp *= calorieToJoule;

//     SpeciesThermoState state;
//     state.volume           = V;
//     state.gibbs_energy     = G;
//     state.enthalpy         = H;
//     state.entropy          = S;
//     state.internal_energy  = U;
//     state.helmholtz_energy = A;
//     state.heat_capacity_cp = Cp;
//     state.heat_capacity_cv = state.heat_capacity_cp; // approximate Cp = Cv for a solid

//     return state;
// }

} // namespace Reaktoro
