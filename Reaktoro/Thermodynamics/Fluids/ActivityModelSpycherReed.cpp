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

#include "ActivityModelSpycherReed.hpp"

// C++ includes
#include <string>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/Utils.hpp>
#include <Reaktoro/Thermodynamics/Solutions/GeneralMixture.hpp>

namespace Reaktoro {

using std::log;

namespace {

// The numbers in the constants and functions below are: 1-H2O, 2-CO2, 3-CH4

// The coefficients for pure gas H2O from Table 1 of Spycher and Reed (1988)
// on the temperature range 0--340 degC and maximum pressure Psat
const auto a11 = -6191.41;
const auto b11 = 14.8528;
const auto c11 = -914.267e-05;
const auto d111 = -6633.26e-02;
const auto e111 = 18277.0e-05;
const auto f111 = -13274.0e-08;

// The coefficients for pure gas CO2 from Table 1 of Spycher and Reed (1988)
// on the temperature range 50--350 C and maximum pressure 500 bar
const auto a22 = -1430.87;
const auto b22 = 3.598;
const auto c22 = -227.376e-05;
const auto d222 = 347.644e-02;
const auto e222 = -1042.47e-05;
const auto f222 = 846.271e-08;

// The coefficients for pure gas CH4 from Table 1 of Spycher and Reed (1988)
// on the temperature range 16--350 C and maximum pressure 500 bar
const auto a33 = -537.779;
const auto b33 = 1.54946;
const auto c33 = -92.7827e-05;
const auto d333 = 120.861e-02;
const auto e333 = -370.814e-05;
const auto f333 = 333.804e-08;

// The coefficients for the binary mixture H2O-CO2 from Table 2 of Spycher and
// Reed (1988) on the temperature range 50--350 C and maximum pressure 94 bar
const auto a12 = -1954.70;
const auto b12 = 7.74805;
const auto c12 = -1.02901e-02;
const auto d112 = 104.453;
const auto e112 = -38.4283e-02;
const auto f112 = 36.5858e-05;
const auto d122 = -8.28426;
const auto e122 = 1.19097e-02;
const auto f122 = 0.808886e-05;

// The coefficients for the binary mixture H2O-CH4 from Table 2 of Spycher and
// Reed (1988) on the temperature range 40--240 C and maximum pressure 500 bar
const auto a13 = -1103.20;
const auto b13 = 4.52871;
const auto c13 = -0.507784e-02;
const auto d113 = 0.0;
const auto e113 = 0.0;
const auto f113 = 0.0;
const auto d133 = 0.0;
const auto e133 = 0.0;
const auto f133 = 0.0;

// The coefficients for the binary mixture CO2-CH4 from Table 2 of Spycher and
// Reed (1988) on the temperature range 25--100 C and maximum pressure 500 bar
const auto a23 = -800.592;
const auto b23 = 2.28990;
const auto c23 = -0.153917e-02;
const auto d223 = 2.99160;
const auto e223 = -1.04893e-02;
const auto f223 = 1.02627e-05;
const auto d233 = 1.58384;
const auto e233 = -0.492077e-02;
const auto f233 = 0.430104e-05;

const auto d123 = 0.0;
const auto e123 = 0.0;
const auto f123 = 0.0;

const double a[][3] =
{
    {a11, a12, a13},
    {a12, a22, a23},
    {a13, a23, a33}
};

const double b[][3] =
{
    {b11, b12, b13},
    {b12, b22, b23},
    {b13, b23, b33}
};

const double c[][3] =
{
    {c11, c12, c13},
    {c12, c22, c23},
    {c13, c23, c33}
};

const double d[][3][3] =
{
    {{d111, d112, d113},
    {d112, d122, d123},
    {d113, d123, d133}},

    {{d112, d122, d123},
    {d122, d222, d223},
    {d123, d223, d233}},

    {{d113, d123, d133},
    {d123, d223, d233},
    {d133, d233, d333}},
};

const double e[][3][3] =
{
    {{e111, e112, e113},
    {e112, e122, e123},
    {e113, e123, e133}},

    {{e112, e122, e123},
    {e122, e222, e223},
    {e123, e223, e233}},

    {{e113, e123, e133},
    {e123, e223, e233},
    {e133, e233, e333}},
};

const double f[][3][3] =
{
    {{f111, f112, f113},
    {f112, f122, f123},
    {f113, f123, f133}},

    {{f112, f122, f123},
    {f122, f222, f223},
    {f123, f223, f233}},

    {{f113, f123, f133},
    {f123, f223, f233},
    {f133, f233, f333}},
};

inline auto computeB(const real& T, int i, int j) -> real
{
    return a[i][j] / (T*T) + b[i][j] / T + c[i][j];
}

inline auto computeBT(const real& T, int i, int j) -> real
{
    return -(2 * a[i][j] / T + b[i][j]) / (T*T);
}

inline auto computeBTT(const real& T, int i, int j) -> real
{
    return (6 * a[i][j] / T + 2 * b[i][j]) / (T*T*T);
}

inline auto computeC(const real& T, int i, int j, int k) -> real
{
    return d[i][j][k] / (T*T) + e[i][j][k] / T + f[i][j][k];
}

inline auto computeCT(const real& T, int i, int j, int k) -> real
{
    return -(2 * d[i][j][k] / T + e[i][j][k]) / (T*T);
}

inline auto computeCTT(const real& T, int i, int j, int k) -> real
{
    return (6 * d[i][j][k] / T + 2 * e[i][j][k]) / (T*T*T);
}

} // namespace

auto fluidChemicalModelSpycherReed(SpeciesListConstRef species)-> ActivityPropsFn
{
    // The names of the gases in the mixture, and the supported ones by this model
    Strings provided = vectorize(species, RKT_LAMBDA(x, x.formula().str()));
    Strings supported = { "H2O", "CO2", "CH4" };

    // Assert the provided gases consists of only a combination of H2O(g), CO2(g) and CH4(g)
    error(!identical(provided, supported), "Could not initialize the "
        "Spycher and Reed (1988) equation of state for the gaseous phase.",
        "This model expects the gaseous phase to be composed of "
        "species H2O(g), CO2(g) and CH4(g) only.");

    // The indices of the species H2O(g), CO2(g) and CH4(g) in the gaseous mixture
    const auto iH2O = species.indexWithFormula("H2O");
    const auto iCO2 = species.indexWithFormula("CO2");
    const auto iCH4 = species.indexWithFormula("CH4");

    // Assert the gaseous species H2O(g), CO2(g) and CH4(g) exist.
    Assert(iH2O < species.size(),
        "Could not create the chemical model Spycher & Reed (1988) for the gaseous phase.",
        "This model requires the species H2O(g) in the gaseous phase.");
    Assert(iCO2 < species.size(),
        "Could not create the chemical model Spycher & Reed (1988) for the gaseous phase.",
        "This model requires the species CO2(g) in the gaseous phase.");
    Assert(iCH4 < species.size(),
        "Could not create the chemical model Spycher & Reed (1988) for the gaseous phase.",
        "This model requires the species CH4(g) in the gaseous phase.");

    // The number of species in the mixture
    const auto nspecies = species.size();

    // The universal gas constant of the phase (in units of J/(mol*K))
    const auto R = universalGasConstant;

    // Define the activity model function of the gaseous phase
    ActivityPropsFn fn = [=](ActivityProps props, ActivityArgs args) mutable
    {
        // The arguments for the activity model evaluation
        const auto& [T, P, x, extra] = args;

        // The pressure in units of bar
        const auto Pbar = 1e-5 * P;

        // The ln of pressure in units of bar
        const auto ln_Pbar = log(Pbar);

        // The ln of mole fractions of the species
        const auto ln_x = x.log();

        // The mole fractions of the gaseous species H2O(g), CO2(g) and CH4(g)
        real y[3] = {};
        if(iH2O < nspecies) y[0] = x[iH2O]; else y[0] = {};
        if(iCO2 < nspecies) y[1] = x[iCO2]; else y[1] = {};
        if(iCH4 < nspecies) y[2] = x[iCH4]; else y[2] = {};

        // Calculate the Bij, BijT, BijTT coefficients
        real B[3][3] = {};
        real BT[3][3] = {};
        real BTT[3][3] = {};
        for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k)
        {
            B[i][k] = computeB(T, i, k);
            BT[i][k] = computeBT(T, i, k);
            BTT[i][k] = computeBTT(T, i, k);
        }

        // Calculate the Cijk, CijkT, CijkTT coefficients
        real C[3][3][3] = {};
        real CT[3][3][3] = {};
        real CTT[3][3][3] = {};
        for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) for (int l = 0; l < 3; ++l)
        {
            C[i][k][l] = computeC(T, i, k, l);
            CT[i][k][l] = computeCT(T, i, k, l);
            CTT[i][k][l] = computeCTT(T, i, k, l);
        }

        // Calculate the coefficient Bmix, BmixT, and BmixTT
        real Bmix = {};
        real BmixT = {};
        real BmixTT = {};
        for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k)
        {
            Bmix += y[i] * y[k] * B[i][k];
            BmixT += y[i] * y[k] * BT[i][k];
            BmixTT += y[i] * y[k] * BTT[i][k];
        }

        // Calculate the coefficient Cmix, CmixT, and CmixTT
        real Cmix = {};
        real CmixT = {};
        real CmixTT = {};
        for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) for (int l = 0; l < 3; ++l)
        {
            Cmix += y[i] * y[k] * y[l] * C[i][k][l];
            CmixT += y[i] * y[k] * y[l] * CT[i][k][l];
            CmixTT += y[i] * y[k] * y[l] * CTT[i][k][l];
        }

        // Calculate the compressibility factor Z and its derivative ZT and ZP
        const real Z = 1.0 + Bmix*Pbar + Cmix*Pbar*Pbar;
        const real ZT = BmixT*Pbar + CmixT*Pbar*Pbar;
        const real ZP = (Bmix + 2*Cmix*Pbar) * 1e-5; // derivative wrt P, not Pbar, that's why 1e-5

        // Calculate the ln fugacity coefficients of the gaseous species
        real ln_phi[3] = {};
        for (int i = 0; i < 3; ++i)
        {
            ln_phi[i] = 1 - Z;
            for (int k = 0; k < 3; ++k)
            {
                ln_phi[i] += 2 * y[k] * B[i][k] * Pbar;
                for (int l = 0; l < 3; ++l)
                    ln_phi[i] += 1.5*y[k] * y[l] * C[i][k][l] * Pbar*Pbar;
            }
        }

        // The result of the chemical model (equation of state) of the phase
        auto& Vex   = props.Vex;
        auto& VexT  = props.VexT;
        auto& VexP  = props.VexP;
        auto& Gres  = props.Gex;
        auto& Hres  = props.Hex;
        auto& Cpres = props.Cpex;
        auto& Cvres = props.Cvex;
        auto& ln_g  = props.ln_g;
        auto& ln_a  = props.ln_a;

        // Calculate the ideal molar volume of the phase (in m3/mol) and its derivatives
        const real V0  =  R*T/P;
        const real V0T =  V0/T;
        const real V0P = -V0/P;

        // Calculate the real volume properties of the phase and its derivatives
        const real V  = Z*V0;
        const real VT = ZT*V0 + Z*V0T;
        const real VP = ZP*V0 + Z*V0P;

        // Calculate the excess molar volume of the phase (in m3/mol) and its derivatives
        Vex  = V - V0;
        VexT = VT - V0T;
        VexP = VP - V0P;

        // Calculate the residual molar Gibbs energy of the phase
        Gres = R * T*(Bmix + 0.5*Cmix*Pbar)*Pbar;

        // Calculate the residual molar enthalpy of the phase
        Hres = -R * T*T*(BmixT + 0.5*CmixT*Pbar)*Pbar;

        // Calculate the residual molar isobaric heat capacity of the phase
        Cpres = 2 * Hres / T - R * T*T*(BmixTT + 0.5*CmixTT*Pbar)*Pbar;

        // Calculate the residual molar isochoric heat capacity of the phase
        Cvres = Cpres + R + T*VT*VT/VP;

        // Set the ln activity coefficients
        ln_g[iH2O] = ln_phi[0];
        ln_g[iCO2] = ln_phi[1];
        ln_g[iCH4] = ln_phi[2];

        // Calculate the ln activities of the species
        ln_a = ln_g + ln_x + ln_Pbar;
    };

    return fn;
}

} // namespace Reaktoro
