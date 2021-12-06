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

#include "AqueousChemicalModelEUNIQUAC.hpp"

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Mixtures/AqueousMixture.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>

namespace Reaktoro {

auto calculateDebyeHuckelParameterA(Temperature T)
{
    const auto coeff0 = 1.131;
    const auto coeff1 = 1.335e-3;
    const auto coeff2 = 1.164e-5;
    return coeff0 + coeff1 * T + coeff2 * T * T;
}

auto aqueousChemicalModelEUNIQUAC(const AqueousMixture& mixture, const EUNIQUACParams& params) -> PhaseChemicalModel
{
    // The molar mass of water
    const double Mw = waterMolarMass;

    // The number of moles of water per kg
    const double nwo = 1/Mw;

    // The number of species in the mixture
    const Index num_species = mixture.numSpecies();

    // The number of charged and neutral species in the mixture
    const Index num_charged_species = mixture.numChargedSpecies();
    const Index num_neutral_species = mixture.numNeutralSpecies();

    // The indices of the charged and neutral species
    const Indices icharged_species = mixture.indicesChargedSpecies();
    const Indices ineutral_species = mixture.indicesNeutralSpecies();

    // The index of the water species
    const Index iwater = mixture.indexWater();

    // The electrical charges of the charged species only
    const Vector charges = mixture.chargesChargedSpecies();

    // The state of the aqueous mixture
    AqueousMixtureState state;

    // Auxiliary variables
    ChemicalScalar xw, ln_xw, sqrtI;

    // Define the intermediate chemical model function of the aqueous mixture
    PhaseChemicalModel model = [=](PhaseChemicalModelResult& res, Temperature T, Pressure P, VectorConstRef n) mutable
    {
        // ******************************************************************************
        // **************** Long-range contribution *************************************
        // ******************************************************************************
        const auto& I = state.Ie;           // ionic strength
        const auto& x = state.x;            // mole fractions of the species
        const auto& m = state.m;            // molalities of the species
        const auto& rho = state.rho/1000;    // density in units of g/cm3

        // Auxiliary references
        auto& ln_g = res.ln_activity_coefficients;
        auto& ln_a = res.ln_activities;

        xw = x[iwater];
        ln_xw = log(xw);
        sqrtI = sqrt(I);
        double b = 1.5;
        auto A_parameter = calculateDebyeHuckelParameterA(T);

        // Loop over all charged species in the mixture
        for(Index i = 0; i < num_charged_species; ++i)
        {
            // The index of the current charged species
            const Index ispecies = icharged_species[i];

            // The molality of the charged species and its molar derivatives
            const auto mi = m[ispecies];

            // The electrical charge of the charged species
            const auto z = charges[i];

            // Calculate the ln activity coefficient of the current charged species
            auto numerator = -z * z * A_parameter * sqrtI;
            auto denominator = 1.0 + b * sqrtI;
            ln_g[ispecies] = numerator / denominator;
        }

        // Computing the water activity coefficient
        auto b_sqrtI = b * sqrtI;
        auto inner_term = 1.0 + b_sqrtI - 1.0 / (1.0 + b_sqrtI) - 2.0 * log(1 + b_sqrtI);
        auto constant_term = Mw * 2.0 * A_parameter / (b * b * b);
        ln_g[iwater] = constant_term * inner_term;

    };

    return model;
}

}
