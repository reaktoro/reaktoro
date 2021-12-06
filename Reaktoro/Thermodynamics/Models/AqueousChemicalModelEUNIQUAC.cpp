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

    // The UNIQUAC parameters ri_values and qi_values of the species
    std::vector<double> ri_values;  // Volume fraction parameter
    std::vector<double> qi_values;  // Surface area parameter

    // Collect the UNIQUAC parameters a and b of the charged species
    for(Index i : icharged_species)
    {
        const AqueousSpecies& species = mixture.species(i);
        ri_values.push_back(params.ri(species.name()));
        qi_values.push_back(params.qi(species.name()));
    }

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
        const double b = 1.5;
        auto A_parameter = calculateDebyeHuckelParameterA(T);

        // Loop over all charged species in the mixture
        for(Index i = 0; i < num_charged_species; ++i)
        {
            // The index of the current charged species
            const Index ispecies = icharged_species[i];

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

        // ******************************************************************************
        // **************** Combinatorial contribution **********************************
        // ******************************************************************************

        // Retrieve water parameters
        const AqueousSpecies& water_species = mixture.species(iwater);
        const auto water_species_name = water_species.name();
        const auto r_w = params.ri(water_species_name);
        const auto q_w = params.qi(water_species_name);

        // Calculate volume fractions (phi_i)

        // First, the auxiliary denominator of phi and theta are computed
        double phi_denominator = 0.0;
        double theta_denominator = 0.0;
        for(Index i = 0; i < num_charged_species; ++i)
        {
            // The index of the current charged species
            const Index ispecies = icharged_species[i];

            // Retrieve UNIQUAC species parameters
            const auto r_i = ri_values[i];
            const auto q_i = qi_values[i];

            // Sum up the charged species i contribution to phi denominator
            const auto xi = x[ispecies];
            phi_denominator += xi.val * r_i;
            theta_denominator += xi.val * q_i;
        }

        // Compute UNIQUAC combinatorial contribution to charged species
        for(Index i = 0; i < num_charged_species; ++i)
        {
            // The index of the current charged species
            const Index ispecies = icharged_species[i];

            // Get species mol fraction
            const auto xi = x[ispecies];

            // Retrieve UNIQUAC species parameters
            const auto r_i = ri_values[i];
            const auto q_i = qi_values[i];

            // Calculate species volume fraction (phi_i)
            auto phi_i = xi * r_i / phi_denominator;

            // Calculate species surface area fraction (theta_i)
            auto theta_i = xi * q_i / theta_denominator;

            // Calculate equation's parts
            auto phi_i_per_xi = phi_i / xi;
            auto ln_phi_i_per_xi = log(phi_i_per_xi);
            auto phi_i_per_theta_i = phi_i / theta_i;
            auto ln_phi_i_per_theta_i = log(phi_i_per_theta_i);

            // Compute the charged species symmetrical combinatorial UNIQUAC contribution
            auto ln_g_combinatorial_sym = ln_phi_i_per_xi + 1.0 - phi_i_per_xi -5.0 * q_i * (ln_phi_i_per_theta_i + 1 - phi_i_per_theta_i);

            // Calculate species combinatorial UNIQUAC activity coeff at infinite dilution.
            // This is necessary to convert the combinatorial contribution to unsymmetrical
            // convention.
            auto ri_rw = r_i / r_w;
            auto qi_qw = q_i / q_w;
            auto ln_g_combinatorial_inf = std::log(ri_rw) + 1.0 - ri_rw;
            ln_g_combinatorial_inf += -5.0 * q_i * (std::log(ri_rw / qi_qw) + 1.0 - ri_rw / qi_qw);

            // Finally, the unsymmetrical combinatorial UNIQUAC contribution
            ln_g[ispecies] += ln_g_combinatorial_sym - ln_g_combinatorial_inf;
        }

    };

    return model;
}

}
