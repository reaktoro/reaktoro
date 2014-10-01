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

#include "AqueousActivityDuanSun.hpp"

// C++ includes
#include <functional>
using namespace std::placeholders;

// Reaktor includes
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Common/ConvertUtils.hpp>

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
        iCO2 = speciesIndex(mixture, "CO2(aq)");
        iNa  = chargedSpeciesLocalIndex(mixture, "Na+");
        iK   = chargedSpeciesLocalIndex(mixture, "K+");
        iCa  = chargedSpeciesLocalIndex(mixture, "Ca++");
        iMg  = chargedSpeciesLocalIndex(mixture, "Mg++");
        iCl  = chargedSpeciesLocalIndex(mixture, "Cl-");
        iSO4 = chargedSpeciesLocalIndex(mixture, "SO4--");
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

auto aqueousActivityDuanSunCO2(const AqueousMixtureState& state, const DuanSunCO2ExtraParams& xparams) -> ChemicalScalar
{
    // Extract temperature and pressure values from the activity parameters
    const double T = state.T;
    const double P = state.P;

    // The molar composition of the aqueous species in the aqueous mixture
    const auto& n = state.n;

    // The molalities of the aqueous species in the aqueous mixture and their molar derivatives
    const auto& m = state.m;

    // The stoichiometric molalities of the ions in the aqueous mixture and their molar derivatives
    const auto& ms = state.ms;

    // The number of species and ions in the aqueous mixture
    const unsigned num_species = n.size();
    const unsigned num_ions = ms.val().size();

    // The parameters lambda and zeta of activity coefficient model
    const double lambda = paramDuanSun(T, P, lambda_coeffs);
    const double zeta   = paramDuanSun(T, P, zeta_coeffs);

    // The zero vector
    const Vector zero = zeros(num_species);

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
    const double mNa_val  = (iNa  < num_ions) ? ms.val().at(iNa)  : 0.0;
    const double mK_val   = (iK   < num_ions) ? ms.val().at(iK)   : 0.0;
    const double mCa_val  = (iCa  < num_ions) ? ms.val().at(iCa)  : 0.0;
    const double mMg_val  = (iMg  < num_ions) ? ms.val().at(iMg)  : 0.0;
    const double mCl_val  = (iCl  < num_ions) ? ms.val().at(iCl)  : 0.0;
    const double mSO4_val = (iSO4 < num_ions) ? ms.val().at(iSO4) : 0.0;

    const Vector mNa_ddn  = (iNa  < num_ions) ? ms.ddn().row(iNa)  : zero;
    const Vector mK_ddn   = (iK   < num_ions) ? ms.ddn().row(iK)   : zero;
    const Vector mCa_ddn  = (iCa  < num_ions) ? ms.ddn().row(iCa)  : zero;
    const Vector mMg_ddn  = (iMg  < num_ions) ? ms.ddn().row(iMg)  : zero;
    const Vector mCl_ddn  = (iCl  < num_ions) ? ms.ddn().row(iCl)  : zero;
    const Vector mSO4_ddn = (iSO4 < num_ions) ? ms.ddn().row(iSO4) : zero;

    // The activity coefficient of CO2(aq) its molar derivatives
    ChemicalScalar gCO2;
    const double gCO2_val = std::exp(2*lambda*(mNa_val + mK_val + 2*mCa_val + 2*mMg_val) +
        zeta*(mNa_val + mK_val + mCa_val + mMg_val)*mCl_val - 0.07*mSO4_val);

    const Vector gCO2_ddn = gCO2_val * (2*lambda*(mNa_ddn + mK_ddn + 2*mCa_ddn + 2*mMg_ddn) +
        zeta*(mNa_ddn + mK_ddn + mCa_ddn + mMg_ddn)*mCl_val +
        zeta*(mNa_val + mK_val + mCa_val + mMg_val)*mCl_ddn - 0.07*mSO4_ddn);

    // The molality of CO2(aq) and its molar derivatives
    const double mCO2_val = m.val().at(iCO2);
    const Vector mCO2_ddn = m.ddn().row(iCO2);

    // The activity of CO2(aq) and its molar derivatives
    const double aCO2_val = mCO2_val * gCO2_val;
    const Vector aCO2_ddn = mCO2_val * gCO2_ddn + mCO2_ddn * gCO2_val;

    return {aCO2_val, 0.0, 0.0, aCO2_ddn};
}

} /* namespace internal */

auto aqueousActivityDuanSunCO2(const AqueousMixture& mixture) -> AqueousActivity
{
    internal::DuanSunCO2ExtraParams xparams(mixture);

    return std::bind(internal::aqueousActivityDuanSunCO2, _1, xparams);
}

} // namespace Reaktor
