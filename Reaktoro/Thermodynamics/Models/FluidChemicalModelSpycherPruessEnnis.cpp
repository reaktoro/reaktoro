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

#include "FluidChemicalModelSpycherPruessEnnis.hpp"

// Reaktoro includes
#include <Reaktoro/Common/ConvertUtils.hpp>
#include <Reaktoro/Math/Roots.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/FluidMixture.hpp>

namespace Reaktoro {
namespace {

// The universal gas constant in units of (bar*cm3)/(mol*K)
const double R = 83.1447;

// Calculates the parameter aCO2 as a function of temperature
inline auto aCO2(Temperature T) -> ThermoScalar
{
    return 7.54e+07 - 4.13e+04 * T;
}

// Parameters from the Table 1 of Spycher et al. (2003)
const double bCO2 = 27.80; // in units of cm3/mol
const double bH2O = 18.18; // in units of cm3/mol
const double aH2OCO2 = 7.89e+07;

/// Calculates the molar volume of the CO2-rich phase (in units of cm3/mol)
auto volumeCO2(Temperature T, ThermoScalar Pb, ThermoScalar sqrtT) -> ThermoScalar
{
    // Auxiliary variables
    const auto amix = aCO2(T);
    const auto bmix = bCO2;

    // The coefficients of the cubic equation
    const ThermoScalar a(1.0);
    const ThermoScalar b = -R * T / Pb;
    const ThermoScalar c = -(R * T * bmix / Pb - amix / (Pb * sqrtT) + bmix * bmix);
    const ThermoScalar d = -amix * bmix / (Pb * sqrtT);

    std::complex<double> x1, x2, x3;
    std::tie(x1, x2, x3) = cardano(a.val, b.val, c.val, d.val);

    double vol = 0;

    if(x2.imag() != 0.0) // there is only one real root
    {
        vol = x1.real();
    } else // there are three real roots
    {
        const auto Vliq = std::min(x1.real(), std::min(x2.real(), x3.real()));
        const auto Vgas = std::max(x1.real(), std::max(x2.real(), x3.real()));

        const auto w1 = Pb * (Vgas - Vliq);
        const auto w2 = R * T * std::log((Vgas - bmix) / (Vliq - bmix)) +
                        amix / (sqrtT * bmix) * std::log((Vgas + bmix) / (Vliq + bmix) * Vliq / Vgas);

        vol = (w2 < w1) ? Vliq : Vgas;
    }

    const double den = 3 * a.val * vol * vol + 2 * b.val * vol + c.val;

    ThermoScalar V;
    V.val = vol;
    V.ddT = -(a.ddT * vol * vol * vol + b.ddT * vol * vol + c.ddT * vol + d.ddT) / den;
    V.ddP = -(a.ddP * vol * vol * vol + b.ddP * vol * vol + c.ddP * vol + d.ddP) / den;

    return V;
}

} // namespace

auto fluidChemicalModelSpycherPruessEnnis(const FluidMixture& mixture) -> PhaseChemicalModel
{
    // The index of the species H2O(g) in the gaseous mixture
    const Index iH2O = mixture.indexSpecies("H2O(g)");

    // The index of the species CO2(g) in the gaseous mixture
    const Index iCO2 = mixture.indexSpecies("CO2(g)");

    // The number of species in the mixture
    const unsigned nspecies = mixture.numSpecies();

    // The ln of H2O(g) and CO2(g) mole fractions
    ChemicalScalar ln_xH2O(nspecies);
    ChemicalScalar ln_xCO2(nspecies);

    // The state of the gaseous mixture
    FluidMixtureState state;

    // Define the chemical model function of the gaseous phase
    PhaseChemicalModel model = [=](PhaseChemicalModelResult& res, Temperature T, Pressure P, VectorConstRef n) mutable {
        // Evaluate the state of the gaseous mixture
        state = mixture.state(T, P, n);

        // Calculate the pressure in bar
        const auto Pb = convertPascalToBar(P);
        const auto ln_Pb = log(Pb);

        // Auxiliary variables
        const auto T05 = sqrt(T);
        const auto T15 = T * T05;

        // Calculate the mixing parameters
        const auto amix = aCO2(T);
        const auto bmix = bCO2;

        // Calculate the molar volume of the CO2-rich phase (in units of cm3/mol)
        const ThermoScalar v = volumeCO2(T, Pb, T05);

        // Auxiliary values for the fugacity coefficients
        const auto aux1 = log(v / (v - bmix));
        const auto aux2 = log((v + bmix) / v) * 2.0 / (R * T15 * bmix);
        const auto aux3 = amix / (R * T15 * bmix * bmix);
        const auto aux4 = log(Pb * v / (R * T));

        // Calculate the fugacity coefficients of H2O(g) and CO2(g) (in natural log scale)
        const auto ln_phiH2O = aux1 + bH2O / (v - bmix) - aH2OCO2 * aux2 +
                               bH2O * aux3 * (log((v + bH2O) / v) - bmix / (v + bmix)) - aux4;

        const auto ln_phiCO2 = aux1 + bCO2 / (v - bmix) - amix * aux2 +
                               bCO2 * aux3 * (log((v + bCO2) / v) - bmix / (v + bmix)) - aux4;

        // The ln mole fractions of all gaseous species
        const ChemicalVector ln_x = log(state.x);

        // The mole fractions of the gaseous species H2O(g) and CO2(g) and their molar derivatives
        if(iH2O < nspecies)
            ln_xH2O = ln_x[iH2O];
        if(iCO2 < nspecies)
            ln_xCO2 = ln_x[iCO2];

        // Set the molar volume of the phase (in units of m3/mol)
        res.molar_volume = convertCubicCentimeterToCubicMeter(v);

        // Set the ln activities of the gaseous species to ideal values
        res.ln_activities = ln_x + ln_Pb;

        // Set the ln activity coefficients of H2O(g) and CO2(g)
        res.ln_activity_coefficients[iH2O] = ln_phiH2O;
        res.ln_activity_coefficients[iCO2] = ln_phiCO2;

        // Correct the ln activities of H2O(g) and CO2(g)
        res.ln_activities[iH2O] += ln_phiH2O;
        res.ln_activities[iCO2] += ln_phiCO2;
    };

    return model;
}

} // namespace Reaktoro
