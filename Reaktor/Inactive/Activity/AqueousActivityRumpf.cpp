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
    : iCO2(mixture.idxSpecies("CO2(aq)")),
      iNa(mixture.idxIon("Na+")),
      iK(mixture.idxIon("K+")),
      iCa(mixture.idxIon("Ca++")),
      iMg(mixture.idxIon("Mg++")),
      iCl(mixture.idxIon("Cl-"))
    {}

    /// The index of the species CO2(aq) in the aqueous mixture
    Index iCO2;

    /// The local index of the ion Na[+] among the ions in a aqueous mixture
    Index iNa;

    /// The local index of the ion K[+] among the ions in a aqueous mixture
    Index iK;

    /// The local index of the ion Ca[2+] among the ions in a aqueous mixture
    Index iCa;

    /// The local index of the ion Mg[2+] among the ions in a aqueous mixture
    Index iMg;

    /// The local index of the ion Cl[-] among the ions in a aqueous mixture
    Index iCl;
};

auto aqueousActivityRumpfCO2(const AqueousActivityParams& params, const RumpfCO2ExtraParams& xparams) -> PartialScalar
{
    // Extract temperature from the parameters
    const double T = params.T;

    // The molar composition of the aqueous solution
    const Vector& n = params.n;

    // The molalities of the aqueous species in the aqueous mixture and their molar derivatives
    const PartialVector& m = params.m;

    // The stoichiometric molalities of the ions in the aqueous mixture and their molar derivatives
    const PartialVector& ms = params.ms;

    // The number of ions in the aqueous mixture
    const double num_ions = func(ms).rows();

    // The index of CO2(aq) in the aqueous mixture
    const Index iCO2 = xparams.iCO2;

    // The local indices of the ions below among all ions in the aqueous mixture
    const Index iNa  = xparams.iNa;
    const Index iK   = xparams.iK;
    const Index iCa  = xparams.iCa;
    const Index iMg  = xparams.iMg;
    const Index iCl  = xparams.iCl;

    // Extract the stoichiometric molalities of the specific ions
    const PartialScalar zero = partialScalar(0.0, zeros(n.rows()));
    const PartialScalar mNa  = (iNa  < num_ions) ? partialScalar(ms, iNa) : zero;
    const PartialScalar mK   = (iK   < num_ions) ? partialScalar(ms, iK)  : zero;
    const PartialScalar mCa  = (iCa  < num_ions) ? partialScalar(ms, iCa) : zero;
    const PartialScalar mMg  = (iMg  < num_ions) ? partialScalar(ms, iMg) : zero;
    const PartialScalar mCl  = (iCl  < num_ions) ? partialScalar(ms, iCl) : zero;

    // The Pitzer's parameters of the Rumpf et al. (1994) model
    const double B = 0.254 - 76.82/T - 10656.0/(T*T) + 6312.0e+3/(T*T*T);
    const double Gamma = -0.0028;

    // The activity coefficient of CO2(aq) its molar derivatives
    PartialScalar gCO2;
    func(gCO2) = std::exp(2*B*(func(mNa) + func(mK) + 2*func(mCa) + 2*func(mMg)) +
        3*Gamma*(func(mNa) + func(mK) + func(mCa) + func(mMg))*func(mCl));

    grad(gCO2) = func(gCO2) * (2*B*(grad(mNa) + grad(mK) + 2*grad(mCa) + 2*grad(mMg)) +
        3*Gamma*(grad(mNa) + grad(mK) + grad(mCa) + grad(mMg))*func(mCl) +
        3*Gamma*(func(mNa) + func(mK) + func(mCa) + func(mMg))*grad(mCl));

    // The molality of CO2(aq) and its molar derivatives
    const PartialScalar mCO2 = partialScalar(m, iCO2);

    // The activity of CO2(aq) and its molar derivatives
    PartialScalar aCO2;
    func(aCO2) = func(mCO2) * func(gCO2);
    grad(aCO2) = func(mCO2) * grad(gCO2) + grad(mCO2) * func(gCO2);

    return aCO2;
}

} /* namespace internal */

auto aqueousActivityRumpfCO2(const AqueousMixture& mixture) -> AqueousActivity
{
    internal::RumpfCO2ExtraParams xparams(mixture);

    return std::bind(internal::aqueousActivityRumpfCO2, _1, xparams);
}

} // namespace Reaktor
