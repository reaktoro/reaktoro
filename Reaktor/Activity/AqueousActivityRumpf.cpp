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
    /// Constructs the instance with provided aqueous solution
    RumpfCO2ExtraParams(const AqueousSolution& solution)
    : iCO2(speciesIndex(solution, "CO2(aq)")),
      iNa(chargedSpeciesLocalIndex(solution, "Na+")),
      iK(chargedSpeciesLocalIndex(solution, "K+")),
      iCa(chargedSpeciesLocalIndex(solution, "Ca++")),
      iMg(chargedSpeciesLocalIndex(solution, "Mg++")),
      iCl(chargedSpeciesLocalIndex(solution, "Cl-"))
    {}

    /// The index of the species CO2(aq) in the aqueous solution
    Index iCO2;

    /// The local index of the ion Na+ among the ions in a aqueous solution
    Index iNa;

    /// The local index of the ion K+ among the ions in a aqueous solution
    Index iK;

    /// The local index of the ion Ca++ among the ions in a aqueous solution
    Index iCa;

    /// The local index of the ion Mg++ among the ions in a aqueous solution
    Index iMg;

    /// The local index of the ion Cl- among the ions in a aqueous solution
    Index iCl;
};

auto computeAqueousActivityRumpfCO2(const AqueousSolutionState& state, const RumpfCO2ExtraParams& xparams) -> ChemicalScalar
{
    // Extract temperature from the parameters
    const double T = state.T;

    // The molar composition of the aqueous solution
    const Vector& n = state.n;

    // The molalities of the aqueous species in the aqueous solution and their molar derivatives
    const ChemicalVector& m = state.m;

    // The stoichiometric molalities of the ions in the aqueous solution and their molar derivatives
    const ChemicalVector& ms = state.ms;

    // The number of species and ions in the aqueous solution
    const unsigned num_species = n.size();
    const unsigned num_ions = ms.val().size();

    // The zero vector
    const Vector zero = zeros(num_species);

    // The index of CO2(aq) in the aqueous solution
    const Index iCO2 = xparams.iCO2;

    // The local indices of the ions below among all ions in the aqueous solution
    const Index iNa  = xparams.iNa;
    const Index iK   = xparams.iK;
    const Index iCa  = xparams.iCa;
    const Index iMg  = xparams.iMg;
    const Index iCl  = xparams.iCl;

    // Extract the stoichiometric molalities of the specific ions and their molar derivatives
    const double mNa_val  = (iNa  < num_ions) ? ms.val()[iNa] : 0.0;
    const double mK_val   = (iK   < num_ions) ? ms.val()[iK]  : 0.0;
    const double mCa_val  = (iCa  < num_ions) ? ms.val()[iCa] : 0.0;
    const double mMg_val  = (iMg  < num_ions) ? ms.val()[iMg] : 0.0;
    const double mCl_val  = (iCl  < num_ions) ? ms.val()[iCl] : 0.0;

    const Vector mNa_ddn  = (iNa  < num_ions) ? ms.ddn().row(iNa) : zero;
    const Vector mK_ddn   = (iK   < num_ions) ? ms.ddn().row(iK)  : zero;
    const Vector mCa_ddn  = (iCa  < num_ions) ? ms.ddn().row(iCa) : zero;
    const Vector mMg_ddn  = (iMg  < num_ions) ? ms.ddn().row(iMg) : zero;
    const Vector mCl_ddn  = (iCl  < num_ions) ? ms.ddn().row(iCl) : zero;

    // The Pitzer's parameters of the Rumpf et al. (1994) model
    const double B = 0.254 - 76.82/T - 10656.0/(T*T) + 6312.0e+3/(T*T*T);
    const double Gamma = -0.0028;

    // The activity coefficient of CO2(aq) its molar derivatives
    const double gCO2_val = std::exp(2*B*(mNa_val + mK_val + 2*mCa_val + 2*mMg_val) +
        3*Gamma*(mNa_val + mK_val + mCa_val + mMg_val)*mCl_val);

    const Vector gCO2_ddn = gCO2_val * (2*B*(mNa_ddn + mK_ddn + 2*mCa_ddn + 2*mMg_ddn) +
        3*Gamma*(mNa_ddn + mK_ddn + mCa_ddn + mMg_ddn)*mCl_val +
        3*Gamma*(mNa_val + mK_val + mCa_val + mMg_val)*mCl_ddn);

    // The molality of CO2(aq) and its molar derivatives
    const double mCO2_val = m.val()[iCO2];
    const Vector mCO2_ddn = m.ddn().row(iCO2);

    // The activity of CO2(aq) and its molar derivatives
    const double aCO2_val = mCO2_val * gCO2_val;
    const Vector aCO2_ddn = mCO2_val * gCO2_ddn + mCO2_ddn * gCO2_val;

    return {aCO2_val, 0.0, 0.0, aCO2_ddn};
}

} // namespace

auto aqueousActivityRumpfCO2(const AqueousSolution& solution) -> AqueousActivity
{
    RumpfCO2ExtraParams xparams(solution);

    return std::bind(computeAqueousActivityRumpfCO2, _1, xparams);
}

} // namespace Reaktor
