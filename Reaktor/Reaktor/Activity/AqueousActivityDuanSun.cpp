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

#include "AqueousActivityDuanSun.hpp"

// C++ includes
#include <functional>
using namespace std::placeholders;

// Reaktor includes
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Common/ConvertUtils.hpp>
#include <Reaktor/Mixtures/AqueousMixture.hpp>

namespace Reaktor {
namespace internal {

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

auto paramDuanSun(double T, double P, const double coeffs[]) -> double
{
    const double Pbar = convert<Pa,bar>(P);

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
        c6*Pbar + c7*Pbar*std::log(T) + c8*Pbar/T + c9*Pbar/(630 - T) +
        c10*Pbar*Pbar/(630 - T)/(630 - T) + c11*T*std::log(Pbar);
}

struct DuanSunCO2ExtraParams
{
    /// Construct the instance with provided aqueous mixture
    DuanSunCO2ExtraParams(const AqueousMixture& mixture)
    {
        iCO2 = indexSpecies(mixture, "CO2(aq)");
        iNa  = localIndexChargedSpecies(mixture, "Na+");
        iK   = localIndexChargedSpecies(mixture, "K+");
        iCa  = localIndexChargedSpecies(mixture, "Ca++");
        iMg  = localIndexChargedSpecies(mixture, "Mg++");
        iCl  = localIndexChargedSpecies(mixture, "Cl-");
        iSO4 = localIndexChargedSpecies(mixture, "SO4--");
    }

    /// The index of the species CO2(aq) in the aqueous mixture
    Index iCO2;

    /// The local index of the ion Na+ among the ions in a aqueous mixture
    Index iNa;

    /// The local index of the ion K+ among the ions in a aqueous mixture
    Index iK;

    /// The local index of the ion Ca++ among the ions in a aqueous mixture
    Index iCa;

    /// The local index of the ion Mg++ among the ions in a aqueous mixture
    Index iMg;

    /// The local index of the ion Cl- among the ions in a aqueous mixture
    Index iCl;

    /// The local index of the ion SO4-- among the ions in a aqueous mixture
    Index iSO4;
};

auto aqueousActivityDuanSunCO2(const AqueousActivityParams& params, const DuanSunCO2ExtraParams& xparams) -> ScalarResult
{
    // Extract temperature and pressure values from the activity parameters
    const double T = params.T;
    const double P = params.P;

    // The molar composition of the aqueous species in the aqueous mixture
    const auto& n = params.n;

    // The molalities of the aqueous species in the aqueous mixture and their molar derivatives
    const auto& m = params.m;

    // The stoichiometric molalities of the ions in the aqueous mixture and their molar derivatives
    const auto& ms = params.ms;

    // The number of ions in the aqueous mixture
    const double num_ions = ms.func.n_rows;

    // The parameters lambda and zeta of activity coefficient model
    const double lambda = paramDuanSun(T, P, lambda_coeffs);
    const double zeta   = paramDuanSun(T, P, zeta_coeffs);

    // The index of CO2(aq) in the aqueous mixture
    const Index iCO2 = xparams.iCO2;

    // The local indices of the ions among all ions in the aqueous mixture
    const Index iNa  = xparams.iNa;
    const Index iK   = xparams.iK;
    const Index iCa  = xparams.iCa;
    const Index iMg  = xparams.iMg;
    const Index iCl  = xparams.iCl;
    const Index iSO4 = xparams.iSO4;

    // The stoichiometric molalities of the specific ions and their molar derivatives
    const ScalarResult zero(0.0, zeros(n.n_rows));
    const ScalarResult mNa  = (iNa  < num_ions) ? ms.row(iNa)  : zero;
    const ScalarResult mK   = (iK   < num_ions) ? ms.row(iK)   : zero;
    const ScalarResult mCa  = (iCa  < num_ions) ? ms.row(iCa)  : zero;
    const ScalarResult mMg  = (iMg  < num_ions) ? ms.row(iMg)  : zero;
    const ScalarResult mCl  = (iCl  < num_ions) ? ms.row(iCl)  : zero;
    const ScalarResult mSO4 = (iSO4 < num_ions) ? ms.row(iSO4) : zero;

    // The activity coefficient of CO2(aq) its molar derivatives
    ScalarResult gCO2;
    gCO2.func = std::exp(2*lambda*(mNa.func + mK.func + 2*mCa.func + 2*mMg.func) +
        zeta*(mNa.func + mK.func + mCa.func + mMg.func)*mCl.func - 0.07*mSO4.func);

    gCO2.grad = gCO2.func * (2*lambda*(mNa.grad + mK.grad + 2*mCa.grad + 2*mMg.grad) +
        zeta*(mNa.grad + mK.grad + mCa.grad + mMg.grad)*mCl.func +
        zeta*(mNa.func + mK.func + mCa.func + mMg.func)*mCl.grad - 0.07*mSO4.grad);

    // The molality of CO2(aq) and its molar derivatives
    ScalarResult mCO2 = m.row(iCO2);

    // The activity of CO2(aq) and its molar derivatives
    ScalarResult aCO2;
    aCO2.func = mCO2.func * gCO2.func;
    aCO2.grad = mCO2.func * gCO2.grad + mCO2.grad * gCO2.func;

    return aCO2;
}

} /* namespace internal */

auto aqueousActivityDuanSunCO2(const AqueousMixture& mixture) -> AqueousActivity
{
    internal::DuanSunCO2ExtraParams xparams(mixture);

    return std::bind(internal::aqueousActivityDuanSunCO2, _1, xparams);
}

} // namespace Reaktor
