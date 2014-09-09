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

#include "AqueousActivityRumpf.hpp"

// C++ includes
#include <functional>
using namespace std::placeholders;

// Reaktor includes
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Mixtures/AqueousMixture.hpp>

namespace Reaktor {
namespace internal {

struct RumpfCO2ExtraParams
{
    /// Constructs the instance with provided aqueous mixture
    RumpfCO2ExtraParams(const AqueousMixture& mixture)
    : iCO2(indexSpecies(mixture, "CO2(aq)")),
      iNa(localIndexChargedSpecies(mixture, "Na+")),
      iK(localIndexChargedSpecies(mixture, "K+")),
      iCa(localIndexChargedSpecies(mixture, "Ca++")),
      iMg(localIndexChargedSpecies(mixture, "Mg++")),
      iCl(localIndexChargedSpecies(mixture, "Cl-"))
    {}

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
};

auto aqueousActivityRumpfCO2(const AqueousActivityParams& params, const RumpfCO2ExtraParams& xparams) -> ScalarResult
{
    // Extract temperature from the parameters
    const double T = params.T;

    // The molar composition of the aqueous solution
    const Vector& n = params.n;

    // The molalities of the aqueous species in the aqueous mixture and their molar derivatives
    const VectorResult& m = params.m;

    // The stoichiometric molalities of the ions in the aqueous mixture and their molar derivatives
    const VectorResult& ms = params.ms;

    // The number of ions in the aqueous mixture
    const double num_ions = ms.func.n_rows;

    // The index of CO2(aq) in the aqueous mixture
    const Index iCO2 = xparams.iCO2;

    // The local indices of the ions below among all ions in the aqueous mixture
    const Index iNa  = xparams.iNa;
    const Index iK   = xparams.iK;
    const Index iCa  = xparams.iCa;
    const Index iMg  = xparams.iMg;
    const Index iCl  = xparams.iCl;

    // Extract the stoichiometric molalities of the specific ions
    const ScalarResult zero(0.0, zeros(n.n_rows));
    const ScalarResult mNa  = (iNa  < num_ions) ? ms.row(iNa) : zero;
    const ScalarResult mK   = (iK   < num_ions) ? ms.row(iK)  : zero;
    const ScalarResult mCa  = (iCa  < num_ions) ? ms.row(iCa) : zero;
    const ScalarResult mMg  = (iMg  < num_ions) ? ms.row(iMg) : zero;
    const ScalarResult mCl  = (iCl  < num_ions) ? ms.row(iCl) : zero;

    // The Pitzer's parameters of the Rumpf et al. (1994) model
    const double B = 0.254 - 76.82/T - 10656.0/(T*T) + 6312.0e+3/(T*T*T);
    const double Gamma = -0.0028;

    // The activity coefficient of CO2(aq) its molar derivatives
    ScalarResult gCO2;
    gCO2.func = std::exp(2*B*(mNa.func + mK.func + 2*mCa.func + 2*mMg.func) +
        3*Gamma*(mNa.func + mK.func + mCa.func + mMg.func)*mCl.func);

    gCO2.grad = gCO2.func * (2*B*(mNa.grad + mK.grad + 2*mCa.grad + 2*mMg.grad) +
        3*Gamma*(mNa.grad + mK.grad + mCa.grad + mMg.grad)*mCl.func +
        3*Gamma*(mNa.func + mK.func + mCa.func + mMg.func)*mCl.grad);

    // The molality of CO2(aq) and its molar derivatives
    const ScalarResult mCO2 = m.row(iCO2);

    // The activity of CO2(aq) and its molar derivatives
    ScalarResult aCO2;
    aCO2.func = mCO2.func * gCO2.func;
    aCO2.grad = mCO2.func * gCO2.grad + mCO2.grad * gCO2.func;

    return aCO2;
}

} /* namespace internal */

auto aqueousActivityRumpfCO2(const AqueousMixture& mixture) -> AqueousActivity
{
    internal::RumpfCO2ExtraParams xparams(mixture);

    return std::bind(internal::aqueousActivityRumpfCO2, _1, xparams);
}

} // namespace Reaktor
