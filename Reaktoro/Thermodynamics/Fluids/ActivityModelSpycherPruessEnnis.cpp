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

#include "ActivityModelSpycherPruessEnnis.hpp"

// Reaktoro includes
#include <Reaktoro/Common/ConvertUtils.hpp>
#include <Reaktoro/Math/Roots.hpp>

namespace Reaktoro {

using std::log;
using std::max;
using std::min;
using std::sqrt;

namespace {

// The universal gas constant in units of (bar*cm3)/(mol*K)
const auto R = 83.1447;

// Calculates the parameter aCO2 as a function of temperature
inline auto aCO2(real T) -> real
{
    return 7.54e+07 - 4.13e+04 * T;
}

// Parameters from the Table 1 of Spycher et al. (2003)
const auto bCO2 = 27.80; // in units of cm3/mol
const auto bH2O = 18.18; // in units of cm3/mol
const auto aH2OCO2 = 7.89e+07;

/// Calculates the molar volume of the CO2-rich phase (in units of cm3/mol)
auto volumeCO2(real T, real Pb, real sqrtT) -> real
{
    // Auxiliary variables
    const auto amix = aCO2(T);
    const auto bmix = bCO2;

    // The coefficients of the cubic equation
    const real b = -R * T / Pb;
    const real c = -(R*T*bmix / Pb - amix / (Pb*sqrtT) + bmix * bmix);
    const real d = -amix * bmix / (Pb*sqrtT);

    const auto [x1, x2, x3] = cardano(b, c, d);

    real vol = {};

    if(x2.imag() != 0.0) // there is only one real root
    {
        vol = x1.real();
    }
    else // there are three real roots
    {
        const auto Vliq = min(x1.real(), min(x2.real(), x3.real()));
        const auto Vgas = max(x1.real(), max(x2.real(), x3.real()));

        const auto w1 = Pb * (Vgas - Vliq);
        const auto w2 = R * T*log((Vgas - bmix) / (Vliq - bmix)) +
            amix / (sqrtT*bmix)*log((Vgas + bmix) / (Vliq + bmix) * Vliq / Vgas);

        vol = (w2 < w1) ? Vliq : Vgas;
    }

    return vol;
}

} // namespace

auto activityModelSpycherPruessEnnis(const SpeciesList& species) -> ActivityModel
{
    // The indices of gases H2O(g) and CO2(g)
    const auto iH2O = species.findWithFormula("H2O");
    const auto iCO2 = species.findWithFormula("CO2");

    // The number of species
    const auto nspecies = species.size();

    // Define the activity model function of the gaseous phase
    ActivityModel fn = [=](ActivityPropsRef props, ActivityArgs args) mutable
    {
        // The arguments for the activity model evaluation
        const auto& [T, P, x] = args;

        // Calculate the pressure in bar
        const auto Pb = convertPascalToBar(P);
        const auto ln_Pb = log(Pb);

        // Check valid ranges for temperature and pressure
        warning(Pb > 600.0, "ActivityModelSpycherPruessEnnis is only valid "
            "for pressures up to 600 bar, but given pressure is ", Pb, " bar!");
        warning(T > 373.15, "ActivityModelSpycherPruessEnnis is only valid for "
            "temperatures up to 100 degC, but given temperature is ", T - 273.15, "degC!");
        warning(T < 285.15, "ActivityModelSpycherPruessEnnis is not valid for "
            "temperatures below 12 degC, but given temperature is ", T - 273.15, "degC!");

        // Auxiliary variables
        const auto T05 = sqrt(T);
        const auto T15 = T * T05;

        // Calculate the mixing parameters
        const auto amix = aCO2(T);
        const auto bmix = bCO2;

        // Calculate the molar volume of the CO2-rich phase (in units of cm3/mol)
        const auto v = volumeCO2(T, Pb, T05);

        // Auxiliary values for the fugacity coefficients
        const auto aux1 = log(v / (v - bmix));
        const auto aux2 = log((v + bmix) / v) * 2.0 / (R*T15*bmix);
        const auto aux3 = amix / (R*T15*bmix*bmix);
        const auto aux4 = log(Pb*v / (R*T));

        // Calculate the fugacity coefficients of H2O(g) and CO2(g) (in natural log scale)
        const auto ln_phiH2O = aux1 + bH2O / (v - bmix) - aH2OCO2 * aux2 +
            bH2O * aux3*(log((v + bH2O) / v) - bmix / (v + bmix)) - aux4;

        const auto ln_phiCO2 = aux1 + bCO2 / (v - bmix) - amix * aux2 +
            bCO2 * aux3*(log((v + bCO2) / v) - bmix / (v + bmix)) - aux4;

        // The ln mole fractions of all gaseous species
        const auto ln_x = x.log();

        // The molar volume of the phase (in m3/mol)
        const auto V = convertCubicCentimeterToCubicMeter(v);

        // Set the corrective molar volume of the phase (in m3/mol)
        props.Vx = V;

        // Set the ln activities of the gaseous species to ideal values
        props.ln_a = ln_x + ln_Pb;

        // Set the ln activity coefficients of H2O(g) and CO2(g)
        if(iH2O < nspecies) props.ln_g[iH2O] = ln_phiH2O;
        if(iCO2 < nspecies) props.ln_g[iCO2] = ln_phiCO2;

        // Correct the ln activities of H2O(g) and CO2(g)
        if(iH2O < nspecies) props.ln_a[iH2O] += ln_phiH2O;
        if(iCO2 < nspecies) props.ln_a[iCO2] += ln_phiCO2;
    };

    return fn;
}

auto ActivityModelSpycherPruessEnnis() -> ActivityModelGenerator
{
    return [](const SpeciesList& species)
    {
        return activityModelSpycherPruessEnnis(species);
    };
}

} // namespace Reaktoro
