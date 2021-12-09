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
    const auto T_ref = 273.15;
    return coeff0 + coeff1 * (T - T_ref) + coeff2 * (T - T_ref) * (T - T_ref);
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

    // Collect the UNIQUAC parameters r_i and q_i of the all species
    for (Index i = 0; i < num_species; ++i)
    {
        const AqueousSpecies& species = mixture.species(i);
        ri_values.push_back(params.ri(species.name()));
        qi_values.push_back(params.qi(species.name()));
    }

    // Build enthalpic BIP matrices for mixture species in the phase
    MatrixXd u_0(num_species, num_species);  // constant enthalpic BIPs term
    MatrixXd u_T(num_species, num_species);  // linear enthalpic BIPs term
    for (Index i = 0; i < num_species; ++i)
    {
        const AqueousSpecies& ispecies = mixture.species(i);
        const auto& ispecies_name = ispecies.name();
        for (Index j = 0; j < num_species; ++j)
        {
            const AqueousSpecies& jspecies = mixture.species(j);
            const auto& jspecies_name = jspecies.name();
            const auto& u_ij_0 = params.uij_0(ispecies_name, jspecies_name);
            const auto& u_ij_T = params.uij_T(ispecies_name, jspecies_name);
            u_0(i, j) = u_ij_0;
            u_T(i, j) = u_ij_T;
        }
    }

    // Auxiliary variables
    ChemicalScalar xw, ln_xw, sqrtI;
    ChemicalVector ln_m;

    // Define the intermediate chemical model function of the aqueous mixture
    PhaseChemicalModel model = [=](PhaseChemicalModelResult& res, Temperature T, Pressure P, VectorConstRef n) mutable
    {
        // Evaluate the state of the aqueous mixture
        auto state = mixture.state(T, P, n);

        // ==============================================================================
        // ================ Long-range contribution =====================================
        // ==============================================================================
        const auto& I = state.Ie;           // ionic strength
        const auto& x = state.x;            // mole fractions of the species
        const auto& m = state.m;            // molalities of the species
        const auto& rho = state.rho/1000;    // density in units of g/cm3

        // Auxiliary references
        auto& ln_g = res.ln_activity_coefficients;
        auto& ln_a = res.ln_activities;

        xw = x[iwater];
        ln_xw = log(xw);
        ln_m = log(m);
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

        // Loop over all neutral species in the mixture
        for(Index i = 0; i < num_neutral_species; ++i)
        {
            // The index of the current neutral species
            const Index ispecies = ineutral_species[i];

            // Calculate the DH ln activity coefficient of the current neutral species
            ln_g[ispecies] = 0.0;
        }

        // Computing the water activity coefficient
        auto b_sqrtI = b * sqrtI;
        auto inner_term = 1.0 + b_sqrtI - 1.0 / (1.0 + b_sqrtI) - 2.0 * log(1 + b_sqrtI);
        auto constant_term = Mw * 2.0 * A_parameter / (b * b * b);
        ln_g[iwater] = constant_term * inner_term;

        // ==============================================================================
        // ================ Combinatorial contribution ==================================
        // ==============================================================================

        // Retrieve water parameters
        const AqueousSpecies& water_species = mixture.species(iwater);
        const auto water_species_name = water_species.name();
        const auto r_w = params.ri(water_species_name);
        const auto q_w = params.qi(water_species_name);

        // First, the auxiliary denominator of phi and theta are computed
        double phi_denominator = 0.0;
        double theta_denominator = 0.0;
        for (Index i = 0; i < num_species; ++i)
        {
            // Retrieve UNIQUAC species parameters
            const auto r_i = ri_values[i];
            const auto q_i = qi_values[i];

            // Sum up the charged species i contribution to phi denominator
            const auto xi = x[i];
            phi_denominator += xi.val * r_i;
            theta_denominator += xi.val * q_i;
        }

        // Compute UNIQUAC combinatorial contribution to all species
        std::vector<double> theta;  // to store theta_i values
        std::vector<double> phi;  // to store phi_i values
        for(Index i = 0; i < num_species; ++i)
        {
            // Get species mol fraction
            const auto xi = x[i];

            // Retrieve UNIQUAC species parameters
            const auto r_i = ri_values[i];
            const auto q_i = qi_values[i];

            // Calculate species volume fraction (phi_i)
            auto phi_i = xi * r_i / phi_denominator;
            phi.push_back(phi_i.val);

            // Calculate species surface area fraction (theta_i)
            auto theta_i = xi * q_i / theta_denominator;
            theta.push_back(theta_i.val);

            // Calculate equation's parts
            auto phi_i_per_xi = phi_i / xi;
            auto ln_phi_i_per_xi = log(phi_i_per_xi);
            auto phi_i_per_theta_i = phi_i / theta_i;
            auto ln_phi_i_per_theta_i = log(phi_i_per_theta_i);

            // Compute the species symmetrical combinatorial UNIQUAC contribution
            auto ln_g_combinatorial_sym = ln_phi_i_per_xi + 1.0 - phi_i_per_xi -5.0 * q_i * (ln_phi_i_per_theta_i + 1.0 - phi_i_per_theta_i);

            // Calculate species combinatorial UNIQUAC activity coeff at infinite dilution.
            // This is necessary to convert the combinatorial contribution to unsymmetrical
            // convention.
            auto ri_rw = r_i / r_w;
            auto qi_qw = q_i / q_w;
            auto ln_g_combinatorial_inf = std::log(ri_rw) + 1.0 - ri_rw;
            ln_g_combinatorial_inf += -5.0 * q_i * (std::log(ri_rw / qi_qw) + 1.0 - ri_rw / qi_qw);

            // Finally, the unsymmetrical combinatorial UNIQUAC contribution
//            ln_g[i] += ln_g_combinatorial_sym - ln_g_combinatorial_inf;
            ln_g[i] += 0.0;
        }

        // ==============================================================================
        // ================ Residual contribution =======================================
        // ==============================================================================

        // Calculate u_ij temperature dependent BIPs. Please note: this is a symmetric matrix.
        MatrixXd u(num_species, num_species);
        for (Index i = 0; i < num_species; ++i)
            for (Index j = 0; j < num_species; ++j)
            {
                const auto& u_ij_0 = u_0(i, j);
                const auto& u_ij_T = u_T(i, j);
                // Note that the temperature T must be given in Kelvin
                const double T_ref = 298.15;
                u(i, j) = u_ij_0 + u_ij_T * (T.val - T_ref);
            }

        // Calculate the enthalpic psi BIPs
        MatrixXd psi(num_species, num_species);
        for (Index i = 0; i < num_species; ++i)
            for (Index j = 0; j < num_species; ++j)
                psi(i, j) = std::exp(-(u(i, j) - u(j, j)) / T.val);

        // Calculate UNIQUAC residual contributions to all species
        for (Index i = 0; i < num_species; ++i)
        {
            // Retrieve UNIQUAC species q_i parameter
            const auto q_i = qi_values[i];

            // Calculate auxiliary summation (theta_l * psi_li product, i index is fixed)
            double theta_psi_product = 0.0;
            for (Index l = 0; l < num_species; ++l)
                theta_psi_product += theta[l] * psi(l, i);

            // Calculate another auxiliary summation (normalized theta_j * psi_ij product)
            double normalized_theta_psi_product = 0.0;
            for (Index j = 0; j < num_species; ++j)
            {
                double denominator = 0.0;
                for (Index l = 0; l < num_species; ++l)
                    denominator += theta[l] * psi(l, j);
                normalized_theta_psi_product += theta[j] * psi(i, j) / denominator;
            }

            // Compute the species symmetrical residual UNIQUAC contribution
            auto ln_g_residual_sym = q_i * (1.0 - std::log(theta_psi_product) - normalized_theta_psi_product);

            // Calculate species residual UNIQUAC activity coeff at infinite dilution.
            // This is necessary to convert the residual contribution to unsymmetrical
            // convention.
            auto ln_g_residual_inf = q_i * (1.0 - std::log(psi(iwater, i)) - psi(i, iwater));

            // Assemble the unsymmetrical residual UNIQUAC contribution
//            ln_g[i] += ln_g_residual_sym - ln_g_residual_inf;
            ln_g[i] += 0.0;
        }

        // ==============================================================================
        // ========== Calculate activities from activities coefficients =================
        // ==============================================================================
        for (Index i = 0; i < num_species; ++i)
            ln_a[i] = ln_g[i] + ln_m[i];

    };

    return model;
}

struct EUNIQUACParams::Impl
{
    /// A map to identify species indices in entalphic BIP matrices
    std::map<std::string, int> bips_species_id_map;

    /// The volume fraction parameters `r_i` of the chemical species.
    std::map<std::string, double> ri_values;

    /// The surface area fraction parameters `q_i` of the chemical species.
    std::map<std::string, double> qi_values;

    /// Matrix to store zeroth order energetic BIP values (Uij_0)
    MatrixXd constant_coeff_bips;

    /// Matrix to store first order energetic BIP values (Uij_t)
    MatrixXd linear_coeff_bips;

    Impl()
    {
    }
};

EUNIQUACParams::EUNIQUACParams()
    : pimpl(new Impl())
{
    setDTUvalues();
}

auto EUNIQUACParams::setDTUvalues() -> void
{
    // Define the species id map to retrieve species indices in BIP matrices
    pimpl->bips_species_id_map = {
        {"H2O(l)", 0},
        {"H+", 1},
        {"Na+", 2},
        {"K+", 3},
        {"NH4+", 4},
        {"Cl-", 5},
        {"SO4--", 6},
        {"HSO4-", 7},
        {"NO3-", 8},
        {"OH-", 9},
        {"CO3--", 10},
        {"HCO3-", 11}
    };

    // Volume fraction parameters
    pimpl->ri_values = {
        {"H2O(l)", 0.92},
        {"H+", 0.13779},
        {"Na+", 1.4034},
        {"K+", 2.2304},
        {"NH4+", 4.8154},
        {"Cl-", 10.386},
        {"SO4--", 12.794},
        {"HSO4-", 19.588},
        {"NO3-", 5.4041},
        {"OH-", 9.3973},
        {"CO3--", 11.09},
        {"HCO3-", 3.755},
        {"NN-", 0.13779},
        {"NN+", 0.13779},
        {"Ca++", 3.87},
        {"Sr++", 4.1852},
        {"Ba++", 1.525}
    };

    // Surface area parameters
    pimpl->qi_values = {
        {"H2O(l)", 1.4},
        {"H+", 1.00E-15},
        {"Na+", 1.199},
        {"K+", 2.4306},
        {"NH4+", 4.6028},
        {"Cl-", 10.197},
        {"SO4--", 12.444},
        {"HSO4-", 22.495},
        {"NO3-", 6.2074},
        {"OH-", 8.8171},
        {"CO3--", 11.32},
        {"HCO3-", 4.712},
        {"NN-", 1.00E-15},
        {"NN+", 1.00E-15},
        {"Ca++", 1.48},
        {"Sr++", 0.2818},
        {"Ba++", 0.2278}
    };

    // Zeroth order energetic BIP values (Uij_0)
    MatrixXd uij_0_DTU {
        {0,      1e5,     733.286, 535.023, 54.0297, 1523.39, 752.879, 602.252, 998.92, 600.495, 328.141,  118.702},
        {1e5,     0,       1.0E+10, 1.0E+10, 1.0E+10, 1.0E+10, 1.0E+10, 1.0E+10, 1.0E+10, 1.0E+10, 1.0E+10, 1.0E+10},
        {733.286, 1.0E+10, 0,       -46.194, 375.977, 1443.23, 845.135, 469.488, 797.474, 1398.14, 476.956, 980.982},
        {535.023, 1.0E+10, -46.194, 0,       184.288, 1465.18, 913.824, 445.673, 818.568, 1805.75, 1370.57, 1123.44},
        {54.0297, 1.0E+10, 375.977, 184.288, 0,       1385.31, 677.178, 418.886, 807.246, 2500,    0,       0},
        {1523.39, 1.0E+10, 1443.23, 1465.18, 1385.31, 2214.81, 2036.06, 0,       2175.02, 1895.52, 2372.94, 2014.18},
        {752.879, 1.0E+10, 845.135, 913.824, 677.178, 2036.06, 1265.83, 1004.16, 1757.79, 1225.67, 1158.91, 956.609},
        {602.252, 1.0E+10, 469.488, 445.673, 418.886, 0,       1004.16, 822.989, 0,       2500,    2500,     2500},
        {998.92,  1.0E+10, 797.474, 818.568, 807.246, 2175.02, 1757.79, 0,       2753.71, 1379.95, 0,       0,},
        {600.495, 1.0E+10, 1398.14, 1805.75, 2500,    1895.52, 1225.67, 2500,    1379.95, 1562.88, 1339.04, 2500},
        {328.141, 1.0E+10, 476.956, 1370.57, 0,       2372.94, 1158.91, 2500,    0,       1339.04, 1065.97, 565.786},
        {118.702, 1.0E+10, 980.982, 1123.44, 0,       2014.18, 956.609, 2500,    0,       2500,    565.786, 253.461},
    };
    pimpl->constant_coeff_bips = uij_0_DTU;

    // First order energetic BIP values (Uij_T)
    MatrixXd uij_T_DTU {
        {0, 0, 0.4872, 0.9936, 0.5855, 14.631, 9.4905, 5.5499, 9.3251, 8.5455, -0.5059, 0.96},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0.4872, 0, 0, 0.119, -0.292, 15.635, 11.681, 6.5329, 9.3047, 20.278, 2.8191, 13.296},
        {0.9936, 0, 0.119, 0, 1.0985, 15.329, 12.278, 5.2231, 10.302, 27.283, 5.0517, 2.7217},
        {0.5855, 0, -0.292, 1.0985, 0, 14.848, 10.356, 5.1053, 9.9705, 0, 0, 0},
        {14.631, 0, 15.635, 15.329, 14.848, 14.436, 12.407, 0, 13.449, 13.628, 8.9402, 11.636},
        {9.4905, 0, 11.681, 12.278, 10.356, 12.407, 8.3194, 8.7877, 6.9013, 8.5902, 6.7255, 8.6656},
        {5.5499, 0, 6.5329, 5.2231, 5.1053, 0, 8.7877, 5.7939, 0, 0, 0, 0},
        {9.3251, 0, 9.3047, 10.302, 9.9705, 13.449, 6.9013, 0, 2.2866, 6.6369, 0, 0},
        {8.5455, 0, 20.278, 27.283, 0, 13.628, 8.5902, 0, 6.6369, 5.6169, 3.8129, 0},
        {-0.5059, 0, 2.8191, 5.0517, 0, 8.9402, 6.7255, 0, 0, 3.8129, -4.4653, -1.2033},
        {0.96, 0, 13.296, 2.7217, 0, 11.636, 8.6656, 0, 0, 0, -1.2033, 1.2824}
    };
    pimpl->linear_coeff_bips = uij_T_DTU;
}

auto EUNIQUACParams::ri(const std::string& name) const -> double
{
    auto it = pimpl->ri_values.find(name);
    // TODO: For now, let's return zero when data is unavailable
    return it != pimpl->ri_values.end() ? it->second : 0.0;
}

auto EUNIQUACParams::qi(const std::string& name) const -> double
{
    auto it = pimpl->qi_values.find(name);
    // TODO: For now, let's return zero when data is unavailable
    return it != pimpl->qi_values.end() ? it->second : 0.0;
}

auto EUNIQUACParams::uij_0(
    const std::string& first_species_name,
    const std::string& second_species_name) const -> double
{
    // TODO: implement exception handling when species names are unavailable
    const auto first_species_id = pimpl->bips_species_id_map.at(first_species_name);
    const auto second_species_id = pimpl->bips_species_id_map.at(second_species_name);
    return pimpl->constant_coeff_bips(first_species_id, second_species_id);
}

auto EUNIQUACParams::uij_T(
    const std::string& first_species_name,
    const std::string& second_species_name) const -> double
{
    // TODO: implement exception handling when species names are unavailable
    const auto first_species_id = pimpl->bips_species_id_map.at(first_species_name);
    const auto second_species_id = pimpl->bips_species_id_map.at(second_species_name);
    return pimpl->linear_coeff_bips(first_species_id, second_species_id);
}
}  // namespace Reaktoro
