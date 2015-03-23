// Reaktor is a C++ library for computational reaction modelling.
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

#include "AqueousActivityRumpf.hpp"

// C++ includes
#include <functional>
using namespace std::placeholders;

// Reaktor includes
#include <Reaktor/Common/Index.hpp>

namespace Reaktor {
namespace {

struct RumpfCO2ExtraParams
{
    /// Constructs the instance with provided aqueous mixture
    RumpfCO2ExtraParams(const AqueousMixture& mixture)
    {
        iCO2 = mixture.indexSpecies("CO2(aq)");
        iNa  = mixture.indexChargedSpecies("Na+");
        iK   = mixture.indexChargedSpecies("K+");
        iCa  = mixture.indexChargedSpecies("Ca++");
        iMg  = mixture.indexChargedSpecies("Mg++");
        iCl  = mixture.indexChargedSpecies("Cl-");
    }

    /// The index of the species CO2(aq) in the aqueous mixture
    Index iCO2;

    /// The local index of the ion Na+ among the ions in the aqueous mixture
    Index iNa;

    /// The local index of the ion K+ among the ions in the aqueous mixture
    Index iK;

    /// The local index of the ion Ca++ among the ions in the aqueous mixture
    Index iCa;

    /// The local index of the ion Mg++ among the ions in the aqueous mixture
    Index iMg;

    /// The local index of the ion Cl- among the ions in the aqueous mixture
    Index iCl;
};

auto computeAqueousActivityRumpfCO2(const AqueousMixtureState& state, const RumpfCO2ExtraParams& xparams) -> ChemicalScalar
{
    // Extract temperature from the parameters
    const double T = state.T;

    // The molar composition of the aqueous mixture
    const Vector& n = state.n;

    // The molalities of the aqueous species in the aqueous mixture and their molar derivatives
    const ChemicalVector& m = state.m;

    // The stoichiometric molalities of the ions in the aqueous mixture and their molar derivatives
    const ChemicalVector& ms = state.ms;

    // The number of species and ions in the aqueous mixture
    const unsigned num_species = n.size();
    const unsigned num_ions = ms.val.size();

    // The index of CO2(aq) in the aqueous mixture
    const Index iCO2 = xparams.iCO2;

    // The local indices of the ions below among all ions in the aqueous mixture
    const Index iNa  = xparams.iNa;
    const Index iK   = xparams.iK;
    const Index iCa  = xparams.iCa;
    const Index iMg  = xparams.iMg;
    const Index iCl  = xparams.iCl;

    // Extract the stoichiometric molalities of the specific ions and their molar derivatives
    ChemicalScalar zero(num_species);
    ChemicalScalar mNa = (iNa < num_ions) ? ms.row(iNa) : zero;
    ChemicalScalar mK  = (iK  < num_ions) ? ms.row(iK)  : zero;
    ChemicalScalar mCa = (iCa < num_ions) ? ms.row(iCa) : zero;
    ChemicalScalar mMg = (iMg < num_ions) ? ms.row(iMg) : zero;
    ChemicalScalar mCl = (iCl < num_ions) ? ms.row(iCl) : zero;

    // The Pitzer's parameters of the Rumpf et al. (1994) model
    const double B = 0.254 - 76.82/T - 10656.0/(T*T) + 6312.0e+3/(T*T*T);
    const double Gamma = -0.0028;

    // The activity coefficient of CO2(aq) its molar derivatives
    ChemicalScalar gCO2;
    gCO2.val = std::exp(2*B*(mNa.val + mK.val + 2*mCa.val + 2*mMg.val) +
        3*Gamma*(mNa.val + mK.val + mCa.val + mMg.val)*mCl.val);

    gCO2.ddn = gCO2.val * (2*B*(mNa.ddn + mK.ddn + 2*mCa.ddn + 2*mMg.ddn) +
        3*Gamma*(mNa.ddn + mK.ddn + mCa.ddn + mMg.ddn)*mCl.val +
        3*Gamma*(mNa.val + mK.val + mCa.val + mMg.val)*mCl.ddn);

    // The molality of CO2(aq) and its molar derivatives
    ChemicalScalar mCO2 = m.row(iCO2);

    // The activity of CO2(aq) and its molar derivatives
    ChemicalScalar aCO2;
    aCO2.val = mCO2.val * gCO2.val;
    aCO2.ddn = mCO2.val * gCO2.ddn + mCO2.ddn * gCO2.val;

    return aCO2;
}

} // namespace

auto aqueousActivityRumpfCO2(const AqueousMixture& mixture) -> AqueousActivityFunction
{
    RumpfCO2ExtraParams xparams(mixture);

    return std::bind(computeAqueousActivityRumpfCO2, _1, xparams);
}

} // namespace Reaktor
