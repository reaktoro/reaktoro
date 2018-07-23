// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

#include <string>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Common/ThermoScalar.hpp>
#include <Reaktoro/Core/Utils.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/GaseousMixture.hpp>
#include <Reaktoro/Thermodynamics/Models/GaseousChemicalModelSpycherReed.hpp>

namespace Reaktoro {
namespace {

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

inline auto computeB(const ThermoScalar& T, int i, int j) -> ThermoScalar
{
    return a[i][j]/(T*T) + b[i][j]/T + c[i][j];
}

inline auto computeBT(const ThermoScalar& T, int i, int j) -> ThermoScalar
{
    return -(2*a[i][j]/T + b[i][j])/(T*T);
}

inline auto computeBTT(const ThermoScalar& T, int i, int j) -> ThermoScalar
{
    return (6*a[i][j]/T + 2*b[i][j])/(T*T*T);
}

inline auto computeC(const ThermoScalar& T, int i, int j, int k) -> ThermoScalar
{
    return d[i][j][k]/(T*T) + e[i][j][k]/T + f[i][j][k];
}

inline auto computeCT(const ThermoScalar& T, int i, int j, int k) -> ThermoScalar
{
    return -(2*d[i][j][k]/T + e[i][j][k])/(T*T);
}

inline auto computeCTT(const ThermoScalar& T, int i, int j, int k) -> ThermoScalar
{
    return (6*d[i][j][k]/T + 2*e[i][j][k])/(T*T*T);
}

} // namespace

auto gaseousChemicalModelSpycherReed(const GaseousMixture& mixture) -> PhaseChemicalModel
{
    // The names of the gases in the mixture, and the supported ones by this model
    std::vector<std::string> provided = names(mixture.species());
    std::vector<std::string> supported = {"H2O(g)", "CO2(g)", "CH4(g)"};

    // Assert the provided gases consists of only a combination of H2O(g), CO2(g) and CH4(g)
    Assert(contained(provided, supported), "Could not initialize the "
        "Spycher and Reed (1988) equation of state for the gaseous phase.",
        "This model expects the gaseous phase to be composed of "
        "species H2O(g), CO2(g) and CH4(g) only.");

    // The indices of the species H2O(g), CO2(g) and CH4(g) in the gaseous mixture
    const Index iH2O = mixture.indexSpecies("H2O(g)");
    const Index iCO2 = mixture.indexSpecies("CO2(g)");
    const Index iCH4 = mixture.indexSpecies("CH4(g)");

    // Assert the gaseous species H2O(g), CO2(g) and CH4(g) exist.
    Assert(iH2O < mixture.numSpecies(),
        "Could not create the chemical model Spycher & Reed (1988) for the gaseous phase.",
        "This model requires the species H2O(g) in the gaseous phase.")
    Assert(iCO2 < mixture.numSpecies(),
        "Could not create the chemical model Spycher & Reed (1988) for the gaseous phase.",
        "This model requires the species CO2(g) in the gaseous phase.")
    Assert(iCH4 < mixture.numSpecies(),
        "Could not create the chemical model Spycher & Reed (1988) for the gaseous phase.",
        "This model requires the species CH4(g) in the gaseous phase.")

    // The number of species in the mixture
    const unsigned nspecies = mixture.numSpecies();

    // An auxiliary zero ChemicalScalar instance
    const ChemicalScalar zero(nspecies);

    // The universal gas constant of the phase (in units of J/(mol*K))
    const double R = universalGasConstant;

    // Define the intermediate chemical model function of the gaseous phase
    auto model = [=](const GaseousMixtureState state)
    {
        // Auxiliary references to state variables
        const auto& T = state.T;
        const auto& P = state.P;
        const auto& x = state.x;

        // The pressure in units of bar
        const auto Pbar = 1e-5 * P;

        // The ln of pressure in units of bar
        const auto ln_Pbar = log(Pbar);

        // The ln of molar fractions of the species
        const auto ln_x = log(x);

        // The molar fractions of the gaseous species H2O(g), CO2(g) and CH4(g)
        ChemicalScalar y[3];
        if(iH2O < nspecies) y[0] = x[iH2O]; else y[0] = zero;
        if(iCO2 < nspecies) y[1] = x[iCO2]; else y[1] = zero;
        if(iCH4 < nspecies) y[2] = x[iCH4]; else y[2] = zero;

        // Calculate the Bij, BijT, BijTT coefficients
        ThermoScalar B[3][3], BT[3][3], BTT[3][3];
        for(int i = 0; i < 3; ++i) for(int k = 0; k < 3; ++k)
        {
            B[i][k] = computeB(T, i, k);
            BT[i][k] = computeBT(T, i, k);
            BTT[i][k] = computeBTT(T, i, k);
        }

        // Calculate the Cijk, CijkT, CijkTT coefficients
        ThermoScalar C[3][3][3], CT[3][3][3], CTT[3][3][3];
        for(int i = 0; i < 3; ++i) for(int k = 0; k < 3; ++k) for(int l = 0; l < 3; ++l)
        {
            C[i][k][l] = computeC(T, i, k, l);
            CT[i][k][l] = computeCT(T, i, k, l);
            CTT[i][k][l] = computeCTT(T, i, k, l);
        }

        // Calculate the coefficient Bmix, BmixT, and BmixTT
        ChemicalScalar Bmix(nspecies), BmixT(nspecies), BmixTT(nspecies);
        for(int i = 0; i < 3; ++i) for(int k = 0; k < 3; ++k)
        {
            Bmix += y[i]*y[k]*B[i][k];
            BmixT += y[i]*y[k]*BT[i][k];
            BmixTT += y[i]*y[k]*BTT[i][k];
        }

        // Calculate the coefficient Cmix, CmixT, and CmixTT
        ChemicalScalar Cmix(nspecies), CmixT(nspecies), CmixTT(nspecies);
        for(int i = 0; i < 3; ++i) for(int k = 0; k < 3; ++k) for(int l = 0; l < 3; ++l)
        {
            Cmix += y[i]*y[k]*y[l]*C[i][k][l];
            CmixT += y[i]*y[k]*y[l]*CT[i][k][l];
            CmixTT += y[i]*y[k]*y[l]*CTT[i][k][l];
        }

        // Calculate the compressibility factor Zmix and its temperature derivative ZmixT
        const ChemicalScalar Zmix = 1.0 + Bmix*Pbar + Cmix*Pbar*Pbar;
        const ChemicalScalar ZmixT = BmixT*Pbar + CmixT*Pbar*Pbar;

        // Calculate the ln fugacity coefficients of the gaseous species
        ChemicalScalar ln_phi[3];
        for(int i = 0; i < 3; ++i)
        {
            ln_phi[i] = 1 - Zmix;
            for(int k = 0; k < 3; ++k)
            {
                ln_phi[i] += 2*y[k]*B[i][k]*Pbar;
                for(int l = 0; l < 3; ++l)
                    ln_phi[i] += 1.5*y[k]*y[l]*C[i][k][l]*Pbar*Pbar;
            }
        }

        // The result of the chemical model (equation of state) of the phase
        PhaseChemicalModelResult res(nspecies);
        auto& V = res.molar_volume;
        auto& GR = res.residual_molar_gibbs_energy;
        auto& HR = res.residual_molar_enthalpy;
        auto& CPR = res.residual_molar_heat_capacity_cp;
        auto& CVR = res.residual_molar_heat_capacity_cv;
        auto& ln_g = res.ln_activity_coefficients;
        auto& ln_c = res.ln_activity_constants;
        auto& ln_a = res.ln_activities;

        // Calculate the molar volume of the phase (in units of m3/mol)
        V = R*T*Zmix/P;

        // Calculate the derivatives dP/dT and dV/dT
        const ChemicalScalar dPdT = P*(1.0/T + ZmixT/Zmix);
        const ChemicalScalar dVdT = V*(1.0/T + ZmixT/Zmix);

        // Calculate the residual molar Gibbs energy of the phase
        GR = R*T*(Bmix + 0.5*Cmix*Pbar)*Pbar;

        // Calculate the residual molar enthalpy of the phase
        HR = -R*T*T*(BmixT + 0.5*CmixT*Pbar)*Pbar;

        // Calculate the residual molar isobaric heat capacity of the phase
        CPR = 2*HR/T - R*T*T*(BmixTT + 0.5*CmixTT*Pbar)*Pbar;

        // Calculate the residual molar isochoric heat capacity of the phase
        CVR = CPR - T*dPdT*dVdT + R;

        // Set the ln activity coefficients
        ln_g[iH2O] = ln_phi[0];
        ln_g[iCO2] = ln_phi[1];
        ln_g[iCH4] = ln_phi[2];

        // Calculate the ln activities of the species
        ln_a = ln_g + ln_x + ln_Pbar;

        // Set the ln activity constants of the species
        ln_c = ln_Pbar;

        return res;
    };

    // Define the chemical model function of the gaseous phase
    PhaseChemicalModel f = [=](double T, double P, const Vector& n)
    {
        // Calculate state of the mixture
        const GaseousMixtureState state = mixture.state(T, P, n);

        return model(state);
    };

    return f;
}

} // namespace Reaktoro
