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

// C++ includes
#include <cmath>
using std::log;

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

auto paramDuanSun(const real& T, const real& P, const real coeffs[]) -> real
{
    // Convert pressure from pascal to bar
    const auto Pbar = 1e-5 * P;

    const auto c1  = coeffs[0];
    const auto c2  = coeffs[1];
    const auto c3  = coeffs[2];
    const auto c4  = coeffs[3];
    const auto c5  = coeffs[4];
    const auto c6  = coeffs[5];
    const auto c7  = coeffs[6];
    const auto c8  = coeffs[7];
    const auto c9  = coeffs[8];
    const auto c10 = coeffs[9];
    const auto c11 = coeffs[10];

    return c1 + c2*T + c3/T + c4*T*T + c5/(630 - T) +
        c6*Pbar + c7*Pbar*log(T) + c8*Pbar/T + c9*Pbar/(630 - T) +
        c10*Pbar*Pbar/(630 - T)/(630 - T) + c11*T*log(Pbar);
}

} // namespace

auto aqueousActivityModelDuanSunCO2(const AqueousMixture& mixture) -> AqueousActivityModel
{
    // The number of speciesn and charged species
    const auto nspecies = mixture.numSpecies();
    const auto nions = mixture.numChargedSpecies();

    // The local indices of some charged species among all charged species
    const auto iNa  = mixture.indexChargedSpeciesAny(alternativeChargedSpeciesNames("Na+"));   // Na+, Na[+]
    const auto iK   = mixture.indexChargedSpeciesAny(alternativeChargedSpeciesNames("K+"));    // K+, K[+]
    const auto iCa  = mixture.indexChargedSpeciesAny(alternativeChargedSpeciesNames("Ca++"));  // Ca++, Ca+2, Ca[+2]
    const auto iMg  = mixture.indexChargedSpeciesAny(alternativeChargedSpeciesNames("Mg++"));  // Mg++, Mg+2, Mg[+2]
    const auto iCl  = mixture.indexChargedSpeciesAny(alternativeChargedSpeciesNames("Cl-"));   // Cl-, Cl[-]
    const auto iSO4 = mixture.indexChargedSpeciesAny(alternativeChargedSpeciesNames("SO4--")); // SO4--, SO4-2, SO4[-2]

    AqueousActivityModel f = [=](const AqueousMixtureState& state) mutable
    {
        const auto T = state.T;
        const auto P = state.P;
        const auto& ms = state.ms;
        const auto lambda = paramDuanSun(T, P, lambda_coeffs);
        const auto zeta   = paramDuanSun(T, P, zeta_coeffs);

        const real mNa  = (iNa  < nions) ? ms[iNa]  : 0.0;
        const real mK   = (iK   < nions) ? ms[iK]   : 0.0;
        const real mCa  = (iCa  < nions) ? ms[iCa]  : 0.0;
        const real mMg  = (iMg  < nions) ? ms[iMg]  : 0.0;
        const real mCl  = (iCl  < nions) ? ms[iCl]  : 0.0;
        const real mSO4 = (iSO4 < nions) ? ms[iSO4] : 0.0;

        // Return the ln activity coefficient of CO2(aq)
        return 2*lambda*(mNa + mK + 2*mCa + 2*mMg) + zeta*(mNa + mK + mCa + mMg)*mCl - 0.07*mSO4;
    };

    return f;
}

} // namespace Reaktoro
