/*
 * Reaktor is a C++ library for computational reaction modelling.
 *
 * Copyright (C) 2014 Allan Leal
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "GaseousActivitySpycherReed.hpp"

// C++ includes
#include <functional>
using namespace std::placeholders;

// Reaktor includes
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Common/ConvertUtils.hpp>
#include <Reaktor/Common/OptimizationUtils.hpp>
#include <Reaktor/Mixtures/GaseousMixture.hpp>

namespace Reaktor {
namespace internal {

// The numbers in the constants and functions below are: 1-H2O, 2-CO2, 3-CH4

// The coefficients for pure gas H2O from Table 1 of Spycher and Reed (1988)
// on the temperature range 0--340 degC and maximum pressure Psat
const double a11  = -6191.41;
const double b11  =  14.8528;
const double c11  = -914.267e-05;
const double d111 = -6633.26e-02;
const double e111 =  18277.0e-05;
const double f111 = -13274.0e-08;

// The coefficients for pure gas CO2 from Table 1 of Spycher and Reed (1988)
// on the temperature range 50--350 C and maximum pressure 500 bar
const double a22  = -1430.87;
const double b22  =  3.598;
const double c22  = -227.376e-05;
const double d222 =  347.644e-02;
const double e222 = -1042.47e-05;
const double f222 =  846.271e-08;

// The coefficients for pure gas CH4 from Table 1 of Spycher and Reed (1988)
// on the temperature range 16--350 C and maximum pressure 500 bar
const double a33  = -537.779;
const double b33  =  1.54946;
const double c33  = -92.7827e-05;
const double d333 =  120.861e-02;
const double e333 = -370.814e-05;
const double f333 =  333.804e-08;

// The coefficients for the binary mixture H2O-CO2 from Table 2 of Spycher and
// Reed (1988) on the temperature range 50--350 C and maximum pressure 94 bar
const double a12  = -1954.70;
const double b12  =  7.74805;
const double c12  = -1.02901e-02;
const double d112 =  104.453;
const double e112 = -38.4283e-02;
const double f112 =  36.5858e-05;
const double d122 = -8.28426;
const double e122 =  1.19097e-02;
const double f122 =  0.808886e-05;

// The coefficients for the binary mixture H2O-CH4 from Table 2 of Spycher and
// Reed (1988) on the temperature range 40--240 C and maximum pressure 500 bar
const double a13  = -1103.20;
const double b13  =  4.52871;
const double c13  = -0.507784e-02;
const double d113 =  0.0;
const double e113 =  0.0;
const double f113 =  0.0;
const double d133 =  0.0;
const double e133 =  0.0;
const double f133 =  0.0;

// The coefficients for the binary mixture CO2-CH4 from Table 2 of Spycher and
// Reed (1988) on the temperature range 25--100 C and maximum pressure 500 bar
const double a23  = -800.592;
const double b23  =  2.28990;
const double c23  = -0.153917e-02;
const double d223 =  2.99160;
const double e223 = -1.04893e-02;
const double f223 =  1.02627e-05;
const double d233 =  1.58384;
const double e233 = -0.492077e-02;
const double f233 =  0.430104e-05;

const double d123 = 0.0;
const double e123 = 0.0;
const double f123 = 0.0;

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

inline auto computeB(double T, int i, int j) -> double
{
    return a[i][j]/(T*T) + b[i][j]/T + c[i][j];
}

inline auto computeC(double T, int i, int j, int k) -> double
{
    return d[i][j][k]/(T*T) + e[i][j][k]/T + f[i][j][k];
}

auto gaseousActivitySpycherReedH2OCO2CH4(const GaseousActivityParams& params, Index iH2O, Index iCO2, Index iCH4) -> std::vector<ThermoScalar>
{
    // The temperature (in units of K) and pressure (in units of bar)
    const double T  = params.T;
    const double Pb = convert<Pa,bar>(params.P);

    // The number of species in the gaseous mixture
    const unsigned num_species = params.n.n_rows;

    // The number of moles of the above gaseous species
    double n[3] = {};
    n[0] = (iH2O < num_species) ? params.n[iH2O] : 0.0;
    n[1] = (iCO2 < num_species) ? params.n[iCO2] : 0.0;
    n[2] = (iCH4 < num_species) ? params.n[iCH4] : 0.0;

    const double nt = n[0] + n[1] + n[2];

    // Calculate the gaseous molar fractions
    const double y[3] = {n[0]/nt, n[1]/nt, n[2]/nt};

    // Calculate the Bij coefficients
    double B[3][3] = {};
    for(int i = 0; i < 3; ++i) for(int k = 0; k < 3; ++k)
        B[i][k] = computeB(T, i, k);

    // Calculate the Cijk coefficients
    double C[3][3][3] = {};
    for(int i = 0; i < 3; ++i) for(int k = 0; k < 3; ++k) for(int l = 0; l < 3; ++l)
        C[i][k][l] = computeC(T, i, k, l);

    // Calculate the molar derivatives of the molar fractions
    double dy[3][3] = {};
    for(int i = 0; i < 3; ++i) for(int k = 0; k < 3; ++k)
        dy[i][k] = ((i == k ? 1.0 : 0.0) - y[i])/nt;

    // Calculate the coefficient Bmix and its derivatives
    double Bmix = 0.0;
    double dBmix[3] = {};

    for(int i = 0; i < 3; ++i) for(int k = 0; k < 3; ++k)
    {
        Bmix += y[i]*y[k]*B[i][k];
        dBmix[0] += (dy[i][0]*y[k] + y[i]*dy[k][0])*B[i][k];
        dBmix[1] += (dy[i][1]*y[k] + y[i]*dy[k][1])*B[i][k];
        dBmix[2] += (dy[i][2]*y[k] + y[i]*dy[k][2])*B[i][k];
    }

    // Calculate the coefficient Cmix and its derivatives
    double Cmix = 0.0;
    double dCmix[3] = {};

    for(int i = 0; i < 3; ++i) for(int k = 0; k < 3; ++k) for(int l = 0; l < 3; ++l)
    {
        Cmix += y[i]*y[k]*y[l]*C[i][k][l];
        dCmix[0] += (dy[i][0]*y[k]*y[l] + y[i]*dy[k][0]*y[l] + y[i]*y[k]*dy[l][0])*C[i][k][l];
        dCmix[1] += (dy[i][1]*y[k]*y[l] + y[i]*dy[k][1]*y[l] + y[i]*y[k]*dy[l][1])*C[i][k][l];
        dCmix[2] += (dy[i][2]*y[k]*y[l] + y[i]*dy[k][2]*y[l] + y[i]*y[k]*dy[l][2])*C[i][k][l];
    }

    // Calculate the fugacity coefficients of the gaseous species
    double phi[3] = {};
    double dphi[3][3] = {};

    for(int i = 0; i < 3; ++i)
    {
        phi[i] -= Bmix*Pb + Cmix*Pb*Pb;
        dphi[i][0] -= dBmix[0]*Pb + dCmix[0]*Pb*Pb;
        dphi[i][1] -= dBmix[1]*Pb + dCmix[1]*Pb*Pb;
        dphi[i][2] -= dBmix[2]*Pb + dCmix[2]*Pb*Pb;

        for(int k = 0; k < 3; ++k)
        {
            phi[i] += 2*y[k]*B[i][k]*Pb;
            dphi[i][0] += 2*dy[k][0]*B[i][k]*Pb;
            dphi[i][1] += 2*dy[k][1]*B[i][k]*Pb;
            dphi[i][2] += 2*dy[k][2]*B[i][k]*Pb;

            for(int l = 0; l < 3; ++l)
            {
                phi[i] += 1.5*y[k]*y[l]*C[i][k][l]*Pb*Pb;
                dphi[i][0] += 1.5*(dy[k][0]*y[l] + y[k]*dy[l][0])*C[i][k][l]*Pb*Pb;
                dphi[i][1] += 1.5*(dy[k][1]*y[l] + y[k]*dy[l][1])*C[i][k][l]*Pb*Pb;
                dphi[i][2] += 1.5*(dy[k][2]*y[l] + y[k]*dy[l][2])*C[i][k][l]*Pb*Pb;
            }
        }

        phi[i] = std::exp(phi[i]);
        dphi[i][0] *= phi[i];
        dphi[i][1] *= phi[i];
        dphi[i][2] *= phi[i];
    }

    // The fugacity coefficients of the gaseous species H2O(g), CO2(g) and CH4(g)
    ThermoScalar phiH2O(phi[0], 0.0, 0.0, zeros(num_species));
    ThermoScalar phiCO2(phi[1], 0.0, 0.0, zeros(num_species));
    ThermoScalar phiCH4(phi[2], 0.0, 0.0, zeros(num_species));

    for(int i = 0; i < 3; ++i)
    {
        if(iH2O < num_species) phiH2O.ddn[iH2O] = dphi[0][0];
        if(iCO2 < num_species) phiH2O.ddn[iCO2] = dphi[0][1];
        if(iCH4 < num_species) phiH2O.ddn[iCH4] = dphi[0][2];

        if(iH2O < num_species) phiCO2.ddn[iH2O] = dphi[1][0];
        if(iCO2 < num_species) phiCO2.ddn[iCO2] = dphi[1][1];
        if(iCH4 < num_species) phiCO2.ddn[iCH4] = dphi[1][2];

        if(iH2O < num_species) phiCH4.ddn[iH2O] = dphi[2][0];
        if(iCO2 < num_species) phiCH4.ddn[iCO2] = dphi[2][1];
        if(iCH4 < num_species) phiCH4.ddn[iCH4] = dphi[2][2];
    }

    // The molar fractions of all gaseous species
    const auto& x = params.x;

    // The molar fractionS of the gaseous species H2O(g), CO2(g)m CH4(g) and their molar derivatives
    ThermoScalar zero = ThermoScalar::zero(num_species);
    ThermoScalar xH2O = (iH2O < num_species) ? x.row(iH2O) : zero;
    ThermoScalar xCO2 = (iCO2 < num_species) ? x.row(iCO2) : zero;
    ThermoScalar xCH4 = (iCH4 < num_species) ? x.row(iCH4) : zero;

    // The activity of the gaseous species H2O(g)
    ThermoScalar aH2O;
    aH2O.val = Pb * (phiH2O.val * xH2O.val);
    aH2O.ddn = Pb * (phiH2O.val * xH2O.ddn + phiH2O.ddn * xH2O.val);

    // The activity of the gaseous species CO2(g)
    ThermoScalar aCO2;
    aCO2.val = Pb * (phiCO2.val * xCO2.val);
    aCO2.ddn = Pb * (phiCO2.val * xCO2.ddn + phiCO2.ddn * xCO2.val);

    // The activity of the gaseous species CH4(g)
    ThermoScalar aCH4;
    aCH4.val = Pb * (phiCH4.val * xCH4.val);
    aCH4.ddn = Pb * (phiCH4.val * xCH4.ddn + phiCH4.ddn * xCH4.val);

    return {aH2O, aCO2, aCH4};
}

} /* namespace internal */

auto gaseousActivitySpycherReedH2OCO2CH4(const GaseousMixture& mixture) -> std::vector<GaseousActivity>
{
    // The index of the species H2O(g) in the gaseous mixture
    const Index iH2O = indexSpecies(mixture, "H2O(g)");

    // The index of the species CO2(g) in the gaseous mixture
    const Index iCO2 = indexSpecies(mixture, "CO2(g)");;

    // The index of the species CH4(g) in the gaseous mixture
    const Index iCH4 = indexSpecies(mixture, "CH4(g)");;

    using functiontype = std::function<decltype(internal::gaseousActivitySpycherReedH2OCO2CH4)>;

    functiontype func(internal::gaseousActivitySpycherReedH2OCO2CH4);

    std::shared_ptr<functiontype> memoized_func = memoizeLastPtr(func);

    std::vector<GaseousActivity> activities(3);
    activities[0] = [=](const GaseousActivityParams& params) { return (*memoized_func)(params, iH2O, iCO2, iCH4)[0]; };
    activities[1] = [=](const GaseousActivityParams& params) { return (*memoized_func)(params, iH2O, iCO2, iCH4)[1]; };
    activities[2] = [=](const GaseousActivityParams& params) { return (*memoized_func)(params, iH2O, iCO2, iCH4)[2]; };

    return activities;
}

} // namespace Reaktor
