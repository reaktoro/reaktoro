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

#include "AqueousActivityModelDuanSunCO2.hpp"

// Reaktoro includes
#include <Reaktoro/Common/NamingUtils.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/AqueousMixture.hpp>

namespace Reaktoro {
namespace {

/// The lambda coefficients for the activity coefficient of CO2(aq)
const double lambda_coeffs[] =
{
    -0.411370585,
     6.07632013e-4,
     97.5347708,
     0.0,
     0.0,
     0.0,
     0.0,
    -0.0237622469,
     0.0170656236,
     0.0,
     1.41335834e-5
};

/// The zeta coefficients for the activity coefficient of CO2(aq)
const double zeta_coeffs[] =
{
     3.36389723e-4,
    -1.98298980e-5,
     0.0,
     0.0,
     0.0,
     0.0,
     0.0,
     2.12220830e-3,
    -5.24873303e-3,
     0.0,
     0.0
};

auto paramDuanSun(const real& T, const real& P, const double coeffs[]) -> real
{
    // Convert pressure from pascal to bar
    const real Pbar = 1e-5 * P;

    const double c1  = coeffs[0];
    const double c2  = coeffs[1];
    const double c3  = coeffs[2];
    const double c4  = coeffs[3];
    const double c5  = coeffs[4];
    const double c6  = coeffs[5];
    const double c7  = coeffs[6];
    const double c8  = coeffs[7];
    const double c9  = coeffs[8];
    const double c10 = coeffs[9];
    const double c11 = coeffs[10];

    return c1 + c2*T + c3/T + c4*T*T + c5/(630 - T) +
        c6*Pbar + c7*Pbar*log(T) + c8*Pbar/T + c9*Pbar/(630 - T) +
        c10*Pbar*Pbar/(630 - T)/(630 - T) + c11*T*log(Pbar);
}

} // namespace

auto aqueousActivityModelDuanSunCO2(const AqueousMixture& mixture) -> AqueousActivityModel
{
    // The number of speciesn and charged species
    const unsigned nspecies = mixture.numSpecies();
    const unsigned nions = mixture.numChargedSpecies();

    // The local indices of some charged species among all charged species
    const Index iNa  = mixture.indexChargedSpeciesAny(alternativeChargedSpeciesNames("Na+"));   // Na+, Na[+]
    const Index iK   = mixture.indexChargedSpeciesAny(alternativeChargedSpeciesNames("K+"));    // K+, K[+]
    const Index iCa  = mixture.indexChargedSpeciesAny(alternativeChargedSpeciesNames("Ca++"));  // Ca++, Ca+2, Ca[+2]
    const Index iMg  = mixture.indexChargedSpeciesAny(alternativeChargedSpeciesNames("Mg++"));  // Mg++, Mg+2, Mg[+2]
    const Index iCl  = mixture.indexChargedSpeciesAny(alternativeChargedSpeciesNames("Cl-"));   // Cl-, Cl[-]
    const Index iSO4 = mixture.indexChargedSpeciesAny(alternativeChargedSpeciesNames("SO4--")); // SO4--, SO4-2, SO4[-2]

    // The molalities of some ionic species covered by the model
    real mNa(nspecies);
    real mK(nspecies);
    real mCa(nspecies);
    real mMg(nspecies);
    real mCl(nspecies);
    real mSO4(nspecies);

    // The ln activity coefficient of CO2(aq)
    real ln_gCO2(nspecies);

    AqueousActivityModel f = [=](const AqueousMixtureState& state) mutable
    {
        // Extract temperature and pressure values
        const real& T = state.T;
        const real& P = state.P;

        // The stoichiometric molalities of the ions in the aqueous mixture and their molar derivatives
        const auto& ms = state.ms;

        // The parameters lambda and zeta of activity coefficient model
        const real lambda = paramDuanSun(T, P, lambda_coeffs);
        const real zeta   = paramDuanSun(T, P, zeta_coeffs);

        // The stoichiometric molalities of the specific ions and their molar derivatives
        if(iNa  < nions) mNa  = ms[iNa];
        if(iK   < nions) mK   = ms[iK];
        if(iCa  < nions) mCa  = ms[iCa];
        if(iMg  < nions) mMg  = ms[iMg];
        if(iCl  < nions) mCl  = ms[iCl];
        if(iSO4 < nions) mSO4 = ms[iSO4];

        // The ln activity coefficient of CO2(aq)
        ln_gCO2 = 2*lambda*(mNa + mK + 2*mCa + 2*mMg) +
            zeta*(mNa + mK + mCa + mMg)*mCl - 0.07*mSO4;

        return ln_gCO2;
    };

    return f;
}

} // namespace Reaktoro
