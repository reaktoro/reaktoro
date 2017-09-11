// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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

#include "MineralChemicalModelIdeal.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/MineralMixture.hpp>

namespace Reaktoro {

auto mineralChemicalModelRedlichKister(const MineralMixture& mixture, double a0, double a1, double a2) -> PhaseChemicalModel
{
    Assert(mixture.numSpecies() == 2,
        "Cannot create the chemical model Redlich-Kister for the mineral phase.",
        "The Redlich-Kister model requires a solid solution phase with two species.");

    // The state of the mineral mixture
    MineralMixtureState state;

    // Define the chemical model function of the mineral phase
    PhaseChemicalModel model = [=](PhaseChemicalModelResult& res, Temperature T, Pressure P, VectorConstRef n) mutable
    {
        // Evaluate the state of the mineral mixture
        state = mixture.state(T, P, n);

        const auto RT = universalGasConstant * state.T;

        const auto x1 = state.x[0];
        const auto x2 = state.x[1];

        res.ln_activity_coefficients[0] = x2*x2*(a0 + a1*(3*x1 - x2) + a2*(x1 - x2)*(5*x1 - x2));
        res.ln_activity_coefficients[1] = x1*x1*(a0 - a1*(3*x2 - x1) + a2*(x2 - x1)*(5*x2 - x1));

        res.ln_activities = res.ln_activity_coefficients + log(state.x);

        res.residual_molar_gibbs_energy = (x1*x2*(a0 + a1*(x1 - x2) + a2*pow((x1 - x2), 2))) * RT;
        res.residual_molar_enthalpy = res.residual_molar_gibbs_energy;
    };

    return model;
}

} // namespace Reaktoro



