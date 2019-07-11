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

#include "SpeciesThermoStateHKF.hpp"

// C++ includes
#include <cmath>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/ConvertUtils.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/GlobalOptions.hpp>
#include <Reaktoro/Common/NamingUtils.hpp>
#include <Reaktoro/Common/ThermoScalar.hpp>
#include <Reaktoro/Thermodynamics/Models/SpeciesElectroState.hpp>
#include <Reaktoro/Thermodynamics/Models/SpeciesElectroStateHKF.hpp>
#include <Reaktoro/Thermodynamics/Models/SpeciesThermoState.hpp>
#include <Reaktoro/Thermodynamics/Species/AqueousSpecies.hpp>
#include <Reaktoro/Thermodynamics/Species/GaseousSpecies.hpp>
#include <Reaktoro/Thermodynamics/Species/FluidSpecies.hpp>
#include <Reaktoro/Thermodynamics/Species/MineralSpecies.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterElectroState.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterElectroStateJohnsonNorton.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterThermoState.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterThermoStateUtils.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterUtils.hpp>

namespace Reaktoro {
namespace {

/// The reference temperature assumed in the HKF equations of state (in units of K)
const double referenceTemperature = 298.15;

/// The reference temperature assumed in the HKF equations of state (in units of bar)
const double referencePressure = 1.0;

/// The reference dielectric constant of water \epsilon
const double referenceDielectricConstant = 78.24385513;

/// The reference Born function Z (dimensionless)
const double referenceBornZ = -1.278055636e-02;

/// The reference Born function Y (dimensionless)
const double referenceBornY = -5.795424563e-05;

/// The reference Born function Q (dimensionless)
const double referenceBornQ =  6.638388994e-12;

/// The reference Born function N (dimensionless)
const double referenceBornN = -2.321814455e-20;

/// The reference Born function U (dimensionless)
const double referenceBornU = 4.872982291e-14;

/// The reference Born function X (dimensionless)
const double referenceBornX = -3.060388224e-07;

/// The \eta constant in the HKF model (in units of (A*cal)/mol)
const double eta = 1.66027e+05;

/// The constant characteristics \Theta of the solvent (in units of K)
const double theta = 228;

/// The constant characteristics \Psi of the solvent (in units of bar)
const double psi = 2600;

template<class SpeciesType>
auto checkTemperatureValidityHKF(Temperature& T, const SpeciesType& species) -> void
{
    // Get the HKF thermodynamic data of the species
    const auto& hkf = species.thermoData().hkf.get();

    // Check if temperature bounds should be enforced
    if(global::options.exception.enforce_temperature_bounds)
        Assert(T >= 0.0 && T <= hkf.Tmax, "Unable to calculate the "
            "thermodynamic properties of species " + species.name() + " using the "
                "revised HKF equations of state.", "The provided temperature `" +
                    std::to_string(T.val) + " K` is either negative or greater than the "
                        "maximum allowed `" + std::to_string(hkf.Tmax) + " K`.");

    // Ensure temperature is not above the maximum allowed
    T = std::min(T.val, hkf.Tmax);
}

auto checkMineralDataHKF(const MineralSpecies& species) -> void
{
    const auto& hkf = species.thermoData().hkf.get();

    const std::string error = "Unable to calculate the thermodynamic properties of mineral species " +
        species.name() + " using the revised HKF equations of state.";

    if(!std::isfinite(hkf.Gf))
        RuntimeError(error, "Missing `Gf` data for this species in the database");

    if(!std::isfinite(hkf.Hf))
        RuntimeError(error, "Missing `Hf` data for this species in the database");

    if(!std::isfinite(hkf.Sr))
        RuntimeError(error, "Missing `Sr` data for this species in the database");

    if(!std::isfinite(hkf.Vr))
        RuntimeError(error, "Missing `Vr` data for this species in the database");
}

} // namespace

auto speciesThermoStateSolventHKF(Temperature T, Pressure P, const WaterThermoState& wt) -> SpeciesThermoState
{
    // Auxiliary data from Helgeson and Kirkham (1974), on page 1098
    const auto Ttr =  273.16;                   // unit: K
    const auto Str =  15.1320 * calorieToJoule; // unit: J/(mol*K)
    const auto Gtr = -56290.0 * calorieToJoule; // unit: J/mol
    const auto Htr = -68767.0 * calorieToJoule; // unit: J/mol
    const auto Utr = -67887.0 * calorieToJoule; // unit: J/mol
    const auto Atr = -55415.0 * calorieToJoule; // unit: J/mol

    const auto Sw = waterMolarMass * wt.entropy;         // unit: J/(mol*K)
    const auto Hw = waterMolarMass * wt.enthalpy;        // unit: J/mol
    const auto Uw = waterMolarMass * wt.internal_energy; // unit: J/mol

    // Calculate the standard molal thermodynamic properties of the aqueous species
    const auto S  = Sw + Str;
    const auto H  = Hw + Htr;
    const auto U  = Uw + Utr;
    const auto G  = Hw - T * (Sw + Str) + Ttr * Str + Gtr;
    const auto A  = Uw - T * (Sw + Str) + Ttr * Str + Atr;
    const auto V  = wt.volume * waterMolarMass;
    const auto Cp = wt.cp * waterMolarMass;
    const auto Cv = wt.cv * waterMolarMass;

    SpeciesThermoState state;
    state.entropy          = S;
    state.enthalpy         = H;
    state.internal_energy  = U;
    state.gibbs_energy     = G;
    state.helmholtz_energy = A;
    state.volume           = V;
    state.heat_capacity_cp = Cp;
    state.heat_capacity_cv = Cv;

    return state;
}

auto speciesThermoStateSoluteHKF(Temperature T, Pressure P, const AqueousSpecies& species, const SpeciesElectroState& aes, const WaterElectroState& wes) -> SpeciesThermoState
{
    // Get the HKF thermodynamic data of the species
    const auto& hkf = species.thermoData().hkf.get();

    // Auxiliary variables
    const auto Pbar = P * 1.0e-05;
    const auto Tr   = referenceTemperature;
    const auto Pr   = referencePressure;
    const auto Zr   = referenceBornZ;
    const auto Yr   = referenceBornY;
    const auto Gf   = hkf.Gf;
    const auto Hf   = hkf.Hf;
    const auto Sr   = hkf.Sr;
    const auto a1   = hkf.a1;
    const auto a2   = hkf.a2;
    const auto a3   = hkf.a3;
    const auto a4   = hkf.a4;
    const auto c1   = hkf.c1;
    const auto c2   = hkf.c2;
    const auto wr   = hkf.wref;
    const auto w    = aes.w;
    const auto wT   = aes.wT;
    const auto wP   = aes.wP;
    const auto wTT  = aes.wTT;
    const auto Z    = wes.bornZ;
    const auto Y    = wes.bornY;
    const auto Q    = wes.bornQ;
    const auto X    = wes.bornX;

    // Calculate the standard molal thermodynamic properties of the aqueous species
    auto V = a1 + a2/(psi + Pbar) +
        (a3 + a4/(psi + Pbar))/(T - theta) - w*Q - (Z + 1)*wP;

    auto G = Gf - Sr*(T - Tr) - c1*(T*log(T/Tr) - T + Tr)
        + a1*(Pbar - Pr) + a2*log((psi + Pbar)/(psi + Pr))
        - c2*((1.0/(T - theta) - 1.0/(Tr - theta))*(theta - T)/theta
        - T/(theta*theta)*log(Tr/T * (T - theta)/(Tr - theta)))
        + 1.0/(T - theta)*(a3*(Pbar - Pr) + a4*log((psi + Pbar)/(psi + Pr)))
        - w*(Z + 1) + wr*(Zr + 1) + wr*Yr*(T - Tr);

    auto H = Hf + c1*(T - Tr) - c2*(1.0/(T - theta) - 1.0/(Tr - theta))
        + a1*(Pbar - Pr) + a2*log((psi + Pbar)/(psi + Pr))
        + (2.0*T - theta)/pow(T - theta, 2)*(a3*(Pbar - Pr)
        + a4*log((psi + Pbar)/(psi + Pr)))
        - w*(Z + 1) + w*T*Y + T*(Z + 1)*wT + wr*(Zr + 1) - wr*Tr*Yr;

    auto S = Sr + c1*log(T/Tr) - c2/theta*(1.0/(T - theta)
        - 1.0/(Tr - theta) + log(Tr/T * (T - theta)/(Tr - theta))/theta)
        + 1.0/pow(T - theta, 2)*(a3*(Pbar - Pr) + a4*log((psi + Pbar)/(psi + Pr)))
        + w*Y + (Z + 1)*wT - wr*Yr;

    auto Cp = c1 + c2/pow(T - theta, 2) - (2.0*T/pow(T - theta, 3))*(a3*(Pbar - Pr)
        + a4*log((psi + Pbar)/(psi + Pr))) + w*T*X + 2.0*T*Y*wT + T*(Z + 1.0)*wTT;

    auto U = H - Pbar*V;

    auto A = U - T*S;

    // Convert the thermodynamic properties of the gas to the standard units
    V  *= calorieToJoule/barToPascal;
    G  *= calorieToJoule;
    H  *= calorieToJoule;
    S  *= calorieToJoule;
    U  *= calorieToJoule;
    A  *= calorieToJoule;
    Cp *= calorieToJoule;

    SpeciesThermoState state;
    state.volume           = V;
    state.gibbs_energy     = G;
    state.enthalpy         = H;
    state.entropy          = S;
    state.internal_energy  = U;
    state.helmholtz_energy = A;
    state.heat_capacity_cp = Cp;
    state.heat_capacity_cv = state.heat_capacity_cp; // approximate Cp = Cv for an aqueous solution

    return state;
}

auto speciesThermoStateHKF(Temperature T, Pressure P, const AqueousSpecies& species) -> SpeciesThermoState
{
    WaterThermoState wt = waterThermoStateWagnerPruss(T, P, StateOfMatter::Liquid);

    if(isAlternativeWaterName(species.name()))
        return speciesThermoStateSolventHKF(T, P, wt);

    WaterElectroState wes = waterElectroStateJohnsonNorton(T, P, wt);

    FunctionG g = functionG(T, P, wt);

    SpeciesElectroState aes = speciesElectroStateHKF(g, species);

    return speciesThermoStateSoluteHKF(T, P, species, aes, wes);
}

auto speciesThermoStateHKF(Temperature T, Pressure P, const FluidSpecies& species) -> SpeciesThermoState
{
    // Check temperature range validity
    checkTemperatureValidityHKF(T, species);

    // Get the HKF thermodynamic data of the species
    const auto& hkf = species.thermoData().hkf.get();

    // Auxiliary variables
    const auto R    = universalGasConstant;
    const auto Pbar = P * 1.0e-5;
    const auto Tr   = referenceTemperature;
    const auto Gf   = hkf.Gf;
    const auto Hf   = hkf.Hf;
    const auto Sr   = hkf.Sr;
    const auto a    = hkf.a;
    const auto b    = hkf.b;
    const auto c    = hkf.c;

    // Calculate the integrals of the heal capacity function of the gas from Tr to T at constant pressure Pr
    const auto CpdT   = a*(T - Tr) + 0.5*b*(T*T - Tr*Tr) - c*(1.0/T - 1.0/Tr);
    const auto CpdlnT = a*log(T/Tr) + b*(T - Tr) - 0.5*c*(1.0/(T*T) - 1.0/(Tr*Tr));

    // Calculate the standard molal thermodynamic properties of the gas
    auto V  = R*T/P; // the ideal gas molar volume (in units of m3/mol)
    auto G  = Gf - Sr * (T - Tr) + CpdT - T * CpdlnT;
    auto H  = Hf + CpdT;
    auto S  = Sr + CpdlnT;
    auto U  = H - Pbar*V;
    auto A  = U - T*S;
    auto Cp = a + b*T + c/(T*T);

    // Convert the thermodynamic properties of the gas to the standard units
    G  *= calorieToJoule;
    H  *= calorieToJoule;
    S  *= calorieToJoule;
    U  *= calorieToJoule;
    A  *= calorieToJoule;
    Cp *= calorieToJoule;

    SpeciesThermoState state;
    state.volume           = V;
    state.gibbs_energy     = G;
    state.enthalpy         = H;
    state.entropy          = S;
    state.internal_energy  = U;
    state.helmholtz_energy = A;
    state.heat_capacity_cp = Cp;
    state.heat_capacity_cv = state.heat_capacity_cp - R;

    return state;
}

auto speciesThermoStateHKF(Temperature T, Pressure P, const MineralSpecies& species) -> SpeciesThermoState
{
    // Check temperature range validity
    checkTemperatureValidityHKF(T, species);

    // Check if the HKF thermodynamic data of the mineral is indeed available
    checkMineralDataHKF(species);

    // Get the HKF thermodynamic data of the species
    const auto& hkf = species.thermoData().hkf.get();

    // Auxiliary variables
    const auto  Pb   = P * 1.0e-5;
    const auto& Tr   = referenceTemperature;
    const auto& Pr   = referencePressure;
    const auto& Gf   = hkf.Gf;
    const auto& Hf   = hkf.Hf;
    const auto& Sr   = hkf.Sr;
    const auto& Vr   = hkf.Vr;
    const auto& nt   = hkf.nptrans;
    const auto& a    = hkf.a;
    const auto& b    = hkf.b;
    const auto& c    = hkf.c;
    const auto& Tt   = hkf.Ttr;
    const auto& dHt  = hkf.Htr;
    const auto& dVt  = hkf.Vtr;
    const auto& dPdT = hkf.dPdTtr;

    // Collect the temperature points used for the integrals along the pressure line P = Pr
    std::vector<Temperature> Ti;

    Ti.push_back(Tr);

    for(int i = 0; i < nt; ++i)
        if(T > Tt[i]) Ti.push_back(Tt[i]);

    Ti.push_back(T);

    // Collect the pressure intercepts along the temperature line T for every phase transition boundary (see
    std::vector<ThermoScalar> Pt;
    for(int i = 0; i < nt; ++i)
    {
        if(dPdT[i] != 0.0)
            Pt.push_back(Pr + dPdT[i]*(T - Tt[i]));
    }

    // Calculate the heat capacity of the mineral at T
    ThermoScalar Cp;
    for(unsigned i = 0; i+1 < Ti.size(); ++i)
        if(Ti[i] <= T && T <= Ti[i+1])
            Cp = a[i] + b[i]*T + c[i]/(T*T);

    // Calculate the integrals of the heat capacity function of the mineral from Tr to T at constant pressure Pr
    ThermoScalar CpdT;
    ThermoScalar CpdlnT;
    for(unsigned i = 0; i+1 < Ti.size(); ++i)
    {
        const auto T0 = Ti[i];
        const auto T1 = Ti[i+1];

        CpdT += a[i]*(T1 - T0) + 0.5*b[i]*(T1*T1 - T0*T0) - c[i]*(1.0/T1 - 1.0/T0);
        CpdlnT += a[i]*log(T1/T0) + b[i]*(T1 - T0) - 0.5*c[i]*(1.0/(T1*T1) - 1.0/(T0*T0));
    }

    // Calculate the volume and other auxiliary quantities for the thermodynamic properties of the mineral
    ThermoScalar V(Vr);
    ThermoScalar GdH;
    ThermoScalar HdH;
    ThermoScalar SdH;
    for(unsigned i = 1; i+1 < Ti.size(); ++i)
    {
        GdH += dHt[i-1]*(T - Ti[i])/Ti[i];
        HdH += dHt[i-1];
        SdH += dHt[i-1]/Ti[i];

        V += dVt[i-1];
    }

    // Calculate the volume integral from Pr to P at constant temperature T
    ThermoScalar VdP = 0.023901488*V*(Pb - Pr);
    for(unsigned i = 0; i < Pt.size(); ++i)
    {
        if(0.0 < Pt[i] && Pt[i] < Pb)
        {
            V   -= dVt[i];
            VdP -= 0.023901488*dVt[i]*(Pb - Pt[i]);
        }
    }

    // Calculate the standard molal thermodynamic properties of the mineral
    auto G = Gf - Sr * (T - Tr) + CpdT - T * CpdlnT + VdP - GdH;
    auto H = Hf + CpdT + VdP + HdH;
    auto S = Sr + CpdlnT + SdH;
    auto U = H - Pb*V;
    auto A = U - T*S;

    // Convert the thermodynamic properties of the mineral to the standard
    V  *= cubicCentimeterToCubicMeter;
    G  *= calorieToJoule;
    H  *= calorieToJoule;
    S  *= calorieToJoule;
    U  *= calorieToJoule;
    A  *= calorieToJoule;
    Cp *= calorieToJoule;

    SpeciesThermoState state;
    state.volume           = V;
    state.gibbs_energy     = G;
    state.enthalpy         = H;
    state.entropy          = S;
    state.internal_energy  = U;
    state.helmholtz_energy = A;
    state.heat_capacity_cp = Cp;
    state.heat_capacity_cv = state.heat_capacity_cp; // approximate Cp = Cv for a solid

    return state;
}

} // namespace Reaktoro
