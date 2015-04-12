// Reaktoro is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "GaseousActivitySpycherPruess.hpp"

// C++ includes
#include <functional>
using namespace std::placeholders;

// Reaktoro includes
#include <Reaktoro/Common/ConvertUtils.hpp>
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/OptimizationUtils.hpp>
#include <Reaktoro/Math/Roots.hpp>

namespace Reaktoro {
namespace {

// The universal gas constant in units of (bar*cm3)/(mol*K)
const double R = 83.1447;

// Calculates the parameter aCO2 as a function of temperature
inline auto aCO2(double T) -> double
{
    return 7.54e+07 - 4.13e+04 * T;
}

// Parameters from the Table 1 of Spycher et al. (2003)
const double bCO2    = 27.80;
const double bH2O    = 18.18;
const double aH2OCO2 = 7.89e+07;

/// Calculates the molar volume of the CO2-rich phase
auto volumeCO2(double T, double Pb, double sqrtT) -> double
{
    // Auxiliary variables
    const double amix = aCO2(T);
    const double bmix = bCO2;

    // The coefficients of the cubic equation
    const double a = 1.0;
    const double b = -R*T/Pb;
    const double c = -(R*T*bmix/Pb - amix/(Pb*sqrtT) + bmix*bmix);
    const double d = -amix*bmix/(Pb*sqrtT);

    std::complex<double> x1, x2, x3;
    std::tie(x1, x2, x3) = cubicRoots(a, b, c, d);

    if(x2.imag() != 0.0) // there is only one real root
    {
        return x1.real();
    }
    else // there are three real roots
    {
        const double Vliq = std::min(x1.real(), std::min(x2.real(), x3.real()));
        const double Vgas = std::max(x1.real(), std::max(x2.real(), x3.real()));

        const double w1 = Pb*(Vgas - Vliq);
        const double w2 = R*T*std::log((Vgas - bmix)/(Vliq - bmix)) +
            amix/(sqrtT*bmix)*std::log((Vgas + bmix)/(Vliq + bmix) * Vliq/Vgas);

        return (w2 < w1) ? Vliq : Vgas;
    }
}

auto computeGaseousActivitiesSpycherPruessH2OCO2(const GaseousMixtureState& state, Index iH2O, Index iCO2) -> std::vector<ChemicalScalar>
{
    // The temperature (in units of K) and pressure (in units of bar)
    const double T  = state.T;
    const double Pb = convert<Pa,bar>(state.P);

    // Auxiliary variables
    const double T05 = std::sqrt(T);
    const double T15 = T * T05;

    const double amix = aCO2(T);
    const double bmix = bCO2;

    // The number of species in the gaseous mixture
    const unsigned num_species = state.n.rows();

    // Calculate the molar volume of the CO2-rich phase
    const double v = volumeCO2(T, Pb, T05);

    // Auxiliary values for the fugacity coefficients
    const double aux1 = std::log(v/(v - bmix));
    const double aux2 = std::log((v + bmix)/v) * 2.0/(R*T15*bmix);
    const double aux3 = amix/(R*T15*bmix*bmix);
    const double aux4 = std::log(Pb*v/(R*T));

    // Calculate the fugacity coefficients of H2O(g) and CO2(g)
    const double phiH2O = std::exp(
        aux1 + bH2O/(v - bmix) - aH2OCO2*aux2 +
        bH2O*aux3*(std::log((v + bH2O)/v) - bmix/(v + bmix)) - aux4);

    const double phiCO2 = std::exp(
        aux1 + bCO2/(v - bmix) - amix*aux2 +
        bCO2*aux3*(std::log((v + bCO2)/v) - bmix/(v + bmix)) - aux4);

    // The molar fractions of all gaseous species
    const auto& x = state.x;

    // The molar fractions of the gaseous species H2O(g) and CO2(g) and their molar derivatives
    ChemicalScalar zero(num_species);
    ChemicalScalar xH2O = (iH2O < num_species) ? x.row(iH2O) : zero;
    ChemicalScalar xCO2 = (iCO2 < num_species) ? x.row(iCO2) : zero;

    // Calculate the activity of the gaseous species H2O(g)
    ChemicalScalar aH2O;
    aH2O.val = phiH2O * Pb * xH2O.val;
    aH2O.ddn = phiH2O * Pb * xH2O.ddn;

    // Calculate the activity of the gaseous species CO2(g)
    ChemicalScalar aCO2;
    aCO2.val = phiCO2 * Pb * xCO2.val;
    aCO2.ddn = phiCO2 * Pb * xCO2.ddn;

    return {aH2O, aCO2};
}

} // namespace

auto gaseousActivitySpycherPruessH2OCO2(const GaseousMixture& mixture) -> std::vector<GaseousActivityFunction>
{
    // The index of the species H2O(g) in the gaseous mixture
    const Index iH2O = mixture.indexSpecies("H2O(g)");

    // The index of the species CO2(g) in the gaseous mixture
    const Index iCO2 = mixture.indexSpecies("CO2(g)");

    using functiontype = std::function<decltype(computeGaseousActivitiesSpycherPruessH2OCO2)>;

    functiontype func(computeGaseousActivitiesSpycherPruessH2OCO2);

    std::shared_ptr<functiontype> memoized_func = memoizeLastPtr(func);

    std::vector<GaseousActivityFunction> activities(2);
    activities[0] = [=](const GaseousMixtureState& params) { return (*memoized_func)(params, iH2O, iCO2)[0]; };
    activities[1] = [=](const GaseousMixtureState& params) { return (*memoized_func)(params, iH2O, iCO2)[1]; };

    return activities;
}

} // namespace Reaktoro
