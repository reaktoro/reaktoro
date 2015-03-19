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

#include "AqueousActivityDrummond.hpp"

// C++ includes
#include <functional>
using namespace std::placeholders;

// Reaktor includes
#include <Reaktor/Common/Index.hpp>

namespace Reaktor {
namespace {

auto computeAqueousActivityDrummondCO2(const AqueousMixtureState& state, Index iCO2) -> ChemicalScalar
{
    // Calculate the activity coefficient of CO2(aq)
    const double T  = state.T;
    const double c1 = -1.0312 + 1.2806e-3*T + 255.9/T;
    const double c2 =  0.4445 - 1.6060e-3*T;

    // The stoichiometric ionic strength of the aqueous mixture
    const auto& I = state.Is;

    // The molalities of the aqueous species
    const auto& m = state.m;

    // The activity coefficient of CO2(aq) and its molar derivatives
    ChemicalScalar gCO2;
    gCO2.val = std::exp(c1 * I.val - c2 * I.val/(I.val + 1));
    gCO2.ddn = gCO2.val * (c1 - c2/((I.val + 1) * (I.val + 1))) * I.ddn;

    // The molality of CO2(aq) and its molar derivatives
    ChemicalScalar mCO2;
    mCO2.val = m.val[iCO2];
    mCO2.ddn = m.ddn.row(iCO2);

    // The activity of CO2(aq) and its molar derivatives
    ChemicalScalar aCO2;
    aCO2.val = mCO2.val * gCO2.val;
    aCO2.ddn = mCO2.val * gCO2.ddn + mCO2.ddn * gCO2.val;

    return aCO2;
}

} // namespace

auto aqueousActivityDrummondCO2(const AqueousMixture& mixture) -> AqueousActivity
{
    Index iCO2 = mixture.indexSpecies("CO2(aq)");

    return std::bind(computeAqueousActivityDrummondCO2, _1, iCO2);
}

} // namespace Reaktor
