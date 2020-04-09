// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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

#include "MineralChemicalModelRedlichKister.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/GeneralMixture.hpp>

namespace Reaktoro {

auto mineralChemicalModelRedlichKister(const GeneralMixture& mixture, double a0, double a1, double a2)-> ActivityPropsFn
{
    error(mixture.numSpecies() != 2, "Cannot create the chemical model Redlich-Kister for the mineral phase. "
        "The Redlich-Kister model requires a solid solution phase with exactly two species.");

    // Define the activity model function of the mineral phase
    ActivityPropsFn fn = [=](ActivityProps res, real T, real P, ArrayXrConstRef x) mutable
    {
        using std::log;
        using std::pow;

        const auto RT = universalGasConstant * T;

        const auto x1 = x[0];
        const auto x2 = x[1];

        res.ln_g[0] = x2*x2*(a0 + a1*(3*x1 - x2) + a2*(x1 - x2)*(5*x1 - x2));
        res.ln_g[1] = x1*x1*(a0 - a1*(3*x2 - x1) + a2*(x2 - x1)*(5*x2 - x1));

        res.ln_a = res.ln_g + log(x);

        res.Gex = (x1*x2*(a0 + a1*(x1 - x2) + a2*pow((x1 - x2), 2))) * RT;
        res.Hex = res.Gex;
    };

    return fn;
}

} // namespace Reaktoro



