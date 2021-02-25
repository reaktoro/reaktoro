// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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

#include "ActivityModelDuanSun.hpp"

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Aqueous/AqueousMixture.hpp>

namespace Reaktoro {
namespace {

using std::log;

/// The lambda coefficients for the activity coefficient of CO2(aq)
const real lambda_coeffs[] =
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
const real zeta_coeffs[] =
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

auto ActivityModelDuanSun(String gas) -> ActivityModel
{
    ActivityModel model = [=](const SpeciesList& species)
    {
        // The index of the dissolved gas in the aqueous phase.
        const auto igas = species.indexWithFormula(gas);

        ActivityPropsFn fn = [=](ActivityPropsRef props, ActivityArgs args)
        {
            // The aqueous mixture and its state exported by a base aqueous activity model.
            const auto& mixture = std::any_cast<AqueousMixture>(args.extra.at(0));
            const auto& state = std::any_cast<AqueousMixtureState>(args.extra.at(1));

            // The local indices of some charged species among all charged species
            static const auto iNa  = mixture.charged().findWithFormula("Na+");
            static const auto iK   = mixture.charged().findWithFormula("K+");
            static const auto iCa  = mixture.charged().findWithFormula("Ca++");
            static const auto iMg  = mixture.charged().findWithFormula("Mg++");
            static const auto iCl  = mixture.charged().findWithFormula("Cl-");
            static const auto iSO4 = mixture.charged().findWithFormula("SO4--");

            const auto& T  = state.T;
            const auto& P  = state.P;
            const auto& ms = state.ms;

            const auto lambda = paramDuanSun(T, P, lambda_coeffs);
            const auto zeta   = paramDuanSun(T, P, zeta_coeffs);

            const auto nions = mixture.charged().size();

            const auto mNa  = (iNa  < nions) ? ms[iNa]  : real(0.0);
            const auto mK   = (iK   < nions) ? ms[iK]   : real(0.0);
            const auto mCa  = (iCa  < nions) ? ms[iCa]  : real(0.0);
            const auto mMg  = (iMg  < nions) ? ms[iMg]  : real(0.0);
            const auto mCl  = (iCl  < nions) ? ms[iCl]  : real(0.0);
            const auto mSO4 = (iSO4 < nions) ? ms[iSO4] : real(0.0);

            props.ln_g[igas] = 2*lambda*(mNa + mK + 2*mCa + 2*mMg) + zeta*(mNa + mK + mCa + mMg)*mCl - 0.07*mSO4;
            props.ln_a[igas] = props.ln_g[igas] + log(state.m[igas]);
        };

        return fn;
    };

    return model;
}

} // namespace Reaktoro
