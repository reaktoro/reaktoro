// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

#include "FluidChemicalModelIdeal.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/FluidMixture.hpp>

namespace Reaktoro {

    auto fluidChemicalModelIdeal(const FluidMixture& mixture) -> PhaseChemicalModel
    {
        // The state of the gaseous mixture
        FluidMixtureState state;

        // Define the chemical model function of the gaseous phase
        PhaseChemicalModel model = [=](PhaseChemicalModelResult& res, Temperature T, Pressure P, VectorConstRef n) mutable
        {
            // Evaluate the state of the gaseous mixture
            state = mixture.state(T, P, n);

            // Calculate pressure in bar
            const ThermoScalar Pbar = 1e-5 * Pressure(P);

            // The ln of pressure in units of bar
            const ThermoScalar ln_Pbar = log(Pbar);

            // The result of the ideal model
            res.ln_activities = log(state.x) + ln_Pbar;
        };

        return model;
    }

} // namespace Reaktoro



