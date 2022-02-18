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
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/AqueousMixture.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>

namespace Reaktoro {

/// Sanity check free function to verify if BIPs matrices have proper dimensions.
/// @see EUNIQUACParams::set_uij_bips
auto sanityCheckEnergyInteractionParams(
    const MatrixXd& uij_0,
    const MatrixXd& uij_T,
    const std::map<std::string, int>& species_id_map) -> void
{
    // Check k's dimensions
    auto uij_0_size = uij_0.size();
    auto uij_T_size = uij_T.size();
    if (uij_0_size > 0 && uij_T_size > 0) {
        auto uij_0_num_of_rows = uij_0.rows();
        auto uij_0_num_of_cols = uij_0.cols();
        auto uij_T_num_of_rows = uij_T.rows();
        auto uij_T_num_of_cols = uij_T.cols();
        auto num_species_in_id_map = species_id_map.size();
        Assert(uij_0_num_of_rows == uij_T_num_of_rows && uij_0_num_of_cols == uij_T_num_of_cols,
            "BIPs matrices should have the same num of rows and columns.", "Check your energy BIPs matrices input."
        );
        Assert(num_species_in_id_map == uij_T_num_of_rows && num_species_in_id_map == uij_T_num_of_cols,
            "The species id map for energy BIPs matrices does not match with provided matrices.",
            "Check your energy BIPs matrices and the related species id map."
        );

        // Check k's symmetry
        for (unsigned i = 0; i < uij_0_num_of_rows; i++){
            for (auto j = i + 1; j < uij_0_num_of_cols; j++){
                Assert(uij_0(i, j) == uij_0(j, i),
                    "BIPs matrix u_ij_0 is not symmetric.", "Check your u_ij_0 BIPs matrix input."
                );
                Assert(uij_T(i, j) == uij_T(j, i),
                    "BIPs matrix u_ij_T is not symmetric.", "Check your u_ij_T BIPs matrix input."
                );
            }
        }
    }
    else {
        Exception exception;
        exception.error << "Invalid dimension for the energy BIPs matrices.";
        exception.reason << "Logic error: a proper non-empty symmetric matrix should be passed as BIPs.";
        RaiseError(exception);
    }
}

auto calculateFittedDebyeHuckelParameterA(Temperature T)
{
    const auto coeff0 = 1.131;
    const auto coeff1 = 1.335e-3;
    const auto coeff2 = 1.164e-5;
    const auto T_ref = 273.15;
    return coeff0 + coeff1 * (T - T_ref) + coeff2 * (T - T_ref) * (T - T_ref);
}

auto calculateDebyeHuckelParameterA(const Temperature T, const ThermoScalar& eps_r, const ThermoScalar& rho)
{
    const auto& F = faradayConstant;
    const auto& R = universalGasConstant;
    const auto& eps_0 = vacuumPermittivity;
    const auto& N_A = avogadroNumber;
    const double pi = 3.14159265359;  // TODO: should we use a more precise value?

    auto constant_factor = F * F * F / (4.0 * pi * N_A);
    auto denominator_factor = eps_0 * eps_r * R * T;
    auto denominator = 2.0 * denominator_factor * denominator_factor * denominator_factor;
    auto sqrt_term = sqrt(rho / denominator);
    auto parameterA = constant_factor * sqrt_term;

    return parameterA;
}

auto aqueousChemicalModelEUNIQUAC(const AqueousMixture& mixture, const EUNIQUACParams& params) -> PhaseChemicalModel
{
    // The molar mass of water
    const double Mw = waterMolarMass;

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
        const auto& rho = state.rho;         // density in units of kg/m3
        const auto& epsilon = state.epsilon; // dielectric constant

        // Auxiliary references
        auto& ln_g = res.ln_activity_coefficients;
        auto& ln_a = res.ln_activities;

        xw = x[iwater];
        ln_xw = log(xw);
        ln_m = log(m);
        sqrtI = sqrt(I);
        const double b = 1.5;
        ThermoScalar A_parameter;
        if (params.useDebyeHuckelGenericParameterA())
            A_parameter = calculateDebyeHuckelParameterA(T, epsilon, rho);
        else
            A_parameter = calculateFittedDebyeHuckelParameterA(T);

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
            ln_g[i] += ln_g_combinatorial_sym - ln_g_combinatorial_inf;
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
            ln_g[i] += ln_g_residual_sym - ln_g_residual_inf;
        }

        // ==============================================================================
        // ========== Calculate activities from activities coefficients =================
        // ==============================================================================

        for (Index i = 0; i < num_charged_species; ++i)
        {
            const Index ispecies = icharged_species[i];
            // Convert activity coefficients to molality scale
            ln_g[ispecies] = ln_g[ispecies] + ln_xw;
            // Activities in molality scale
            ln_a[ispecies] = ln_g[ispecies] + ln_m[ispecies];
        }

        for (Index i = 0; i < num_neutral_species; ++i)
        {
            const Index ispecies = ineutral_species[i];
            // Convert activity coefficients to molality scale
            ln_g[ispecies] = ln_g[ispecies] + ln_xw;
            // Activities in molality scale
            ln_a[ispecies] = ln_g[ispecies] + ln_m[ispecies];
        }

        // Finalize the computation of the activity of water (in mole fraction scale)
        ln_a[iwater] = ln_g[iwater] + ln_xw;

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

    /// Set if Debye-Huckel solvent A-parameter is the fitted expression or the general is used instead.
    bool useGeneralDebyeHuckelParameterA;

    Impl()
    {
    }
};

EUNIQUACParams::EUNIQUACParams()
    : pimpl(new Impl())
{
    setDTUvalues();
    pimpl->useGeneralDebyeHuckelParameterA = false;
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

auto EUNIQUACParams::setVillafafilaGarcia2006() -> void
{
    // Define the species id map to retrieve species indices in BIP matrices
    pimpl->bips_species_id_map = {
        {"H2O(l)", 0},
        {"CO2(aq)", 1},
        {"Na+", 2},
        {"H+", 3},
        {"Ba++", 4},
        {"Ca++", 5},
        {"Sr++", 6},
        {"Mg++", 7},
        {"OH-", 8},
        {"Cl-", 9},
        {"SO4--", 10},
        {"CO3--", 11},
        {"HCO3-", 12}
    };

    // Volume fraction parameters
    pimpl->ri_values = {
        {"H2O(l)", 0.92},
        {"CO2(aq)", 0.75},
        {"Na+", 1.40},
        {"H+", 0.14},
        {"Ba++", 15.67},
        {"Ca++", 3.87},
        {"Sr++", 7.14},
        {"Mg++", 5.41},
        {"OH-", 9.40},
        {"Cl-", 10.39},
        {"SO4--", 12.79},
        {"CO3--", 10.83},
        {"HCO3-", 8.08}
    };

    // Surface area parameters
    pimpl->qi_values = {
        {"H2O(l)", 1.4},
        {"CO2(aq)", 2.45},
        {"Na+", 1.2},
        {"H+", 1e-16},
        {"Ba++", 14.48},
        {"Ca++", 1.48},
        {"Sr++", 12.89},
        {"Mg++", 2.54},
        {"OH-", 8.88},
        {"Cl-", 10.20},
        {"SO4--", 12.44},
        {"CO3--", 10.77},
        {"HCO3-", 8.68}
    };

    // Initialize the BIPs as zero matrices to be filled
    const auto& id_map = pimpl->bips_species_id_map;
    auto num_species = pimpl->bips_species_id_map.size();
    MatrixXd uij_0;
    MatrixXd uij_T;
    uij_0.setZero(num_species, num_species);
    uij_T.setZero(num_species, num_species);

    // Auxiliary variables
    std::string species_name_1, species_name_2;
    double bip_value;

    // Auxiliary function to define BIPs symmetrically
    auto set_bip = [] (
        MatrixXd& bips,
        const std::map<std::string, int>& id_map,
        const std::string& species_1,
        const std::string& species_2,
        double value) -> void {
        auto id_species_1 = id_map.at(species_1);
        auto id_species_2 = id_map.at(species_2);
        bips(id_species_1, id_species_2) = value;
        bips(id_species_2, id_species_1) = value;
    };

    // ====================================================================
    // Zeroth order energetic BIP values (Uij_0)
    // ====================================================================
    // Water
    species_name_1 = "H2O(l)";
    set_bip(uij_0, id_map, species_name_1, species_name_1, 0.0);

    species_name_2 = "CO2(aq)";
    bip_value = 8.83825;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "H+";
    bip_value = 1e4;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Na+";
    bip_value = 733.2863;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Ca++";
    bip_value = 496.3523;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Ba++";
    bip_value = -0.3786;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Sr++";
    bip_value = 543.11;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Mg++";
    bip_value = -2.04282;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "OH-";
    bip_value = 600.4952;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Cl-";
    bip_value = 1523.393;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "SO4--";
    bip_value = 752.8792;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "CO3--";
    bip_value = 361.3877;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "HCO3-";
    bip_value = 577.0502;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    // CO2

    species_name_1 = "CO2(aq)";
    set_bip(uij_0, id_map, species_name_1, species_name_1, 302.25);

    species_name_2 = "H+";
    bip_value = 1e10;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Na+";
    bip_value = 172.39;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Ca++";
    bip_value = 2500.0;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Ba++";
    bip_value = 2500.0;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Sr++";
    bip_value = -100.7;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Mg++";
    bip_value = -581.2;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "OH-";
    bip_value = 2500.0;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Cl-";
    bip_value = 1613.0;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "SO4--";
    bip_value = 1942.4;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "CO3--";
    bip_value = 2500.0;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "HCO3-";
    bip_value = 526.31;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    // H+

    species_name_1 = "H+";
    set_bip(uij_0, id_map, species_name_1, species_name_1, 0.0);

    species_name_2 = "Na+";
    bip_value = 1e10;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Ca++";
    bip_value = 1e10;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Ba++";
    bip_value = 1e10;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Sr++";
    bip_value = 1e10;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Mg++";
    bip_value = 1e10;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "OH-";
    bip_value = 1e10;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Cl-";
    bip_value = 1e10;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "SO4--";
    bip_value = 1e10;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "CO3--";
    bip_value = 1e10;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "HCO3-";
    bip_value = 1e10;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    // Na+

    species_name_1 = "Na+";
    set_bip(uij_0, id_map, species_name_1, species_name_1, 0.0);

    species_name_2 = "Ca++";
    bip_value = -100.0;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Ba++";
    bip_value = 779.1;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Sr++";
    bip_value = -103.9;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Mg++";
    bip_value = -70.956;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "OH-";
    bip_value = 1398.1;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Cl-";
    bip_value = 1443.2;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "SO4--";
    bip_value = 845.14;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "CO3--";
    bip_value = 547.95;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "HCO3-";
    bip_value = 1101.9;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    // Ca++

    species_name_1 = "Ca++";
    set_bip(uij_0, id_map, species_name_1, species_name_1, 0.0);

    species_name_2 = "Ba++";
    bip_value = 2500.0;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Sr++";
    bip_value = -402.78;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Mg++";
    bip_value = 155.2324;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "OH-";
    bip_value = 164.6378;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Cl-";
    bip_value = 1805.59;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "SO4--";
    bip_value = 1258.103;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "CO3--";
    bip_value = 2500.0;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "HCO3-";
    bip_value = 2500.0;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    // Ba++

    species_name_1 = "Ba++";
    set_bip(uij_0, id_map, species_name_1, species_name_1, 0.0);

    species_name_2 = "Sr++";
    bip_value = 2500.0;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Mg++";
    bip_value = 628.529;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "OH-";
    bip_value = 2500.0;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Cl-";
    bip_value = 1403.17;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "SO4--";
    bip_value = 2500.0;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "CO3--";
    bip_value = 2500.0;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "HCO3-";
    bip_value = 2500.0;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    // Sr++

    species_name_1 = "Sr++";
    set_bip(uij_0, id_map, species_name_1, species_name_1, 0.0);

    species_name_2 = "Mg++";
    bip_value = -400.58;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "OH-";
    bip_value = 2500.0;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Cl-";
    bip_value = 1895.88;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "SO4--";
    bip_value = 2500.0;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "CO3--";
    bip_value = 2500.0;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "HCO3-";
    bip_value = 2500.0;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    // Mg++

    species_name_1 = "Mg++";
    set_bip(uij_0, id_map, species_name_1, species_name_1, 0.0);

    species_name_2 = "OH-";
    bip_value = 736.42;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Cl-";
    bip_value = 2049.0;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "SO4--";
    bip_value = 1407.21;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "CO3--";
    bip_value = 100.0;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "HCO3-";
    bip_value = 100.0;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    // OH-

    species_name_1 = "OH-";
    set_bip(uij_0, id_map, species_name_1, species_name_1, 1562.88);

    species_name_2 = "Cl-";
    bip_value = 1895.52;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "SO4--";
    bip_value = 1225.67;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "CO3--";
    bip_value = 1588.03;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "HCO3-";
    bip_value = 2500.0;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    // Cl-

    species_name_1 = "Cl-";
    set_bip(uij_0, id_map, species_name_1, species_name_1, 2214.81);

    species_name_2 = "SO4--";
    bip_value = 2036.06;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "CO3--";
    bip_value = 2724.94;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "HCO3-";
    bip_value = 1736.62;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    // SO4--

    species_name_1 = "SO4--";
    set_bip(uij_0, id_map, species_name_1, species_name_1, 1265.83);

    species_name_2 = "CO3--";
    bip_value = 1216.76;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "HCO3-";
    bip_value = 990.48;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    // CO3--

    species_name_1 = "CO3--";
    set_bip(uij_0, id_map, species_name_1, species_name_1, 1458.344);

    species_name_2 = "HCO3-";
    bip_value = 800.0081;
    set_bip(uij_0, id_map, species_name_1, species_name_2, bip_value);

    // HCO3-

    species_name_1 = "HCO3-";
    set_bip(uij_0, id_map, species_name_1, species_name_1, 771.038);

    pimpl->constant_coeff_bips = uij_0;

    // ====================================================================
    // First order energetic BIP values (Uij_T)
    // ====================================================================
    // Water
    species_name_1 = "H2O(l)";
    set_bip(uij_T, id_map, species_name_1, species_name_1, 0.0);

    species_name_2 = "CO2(aq)";
    bip_value = 0.86293;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "H+";
    bip_value = 0.0;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Na+";
    bip_value = 0.48719;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Ca++";
    bip_value = -8.0654;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Ba++";
    bip_value = 0.58244;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Sr++";
    bip_value = 1.2742;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Mg++";
    bip_value = -3.5542;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "OH-";
    bip_value = 8.5455;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Cl-";
    bip_value = 14.631;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "SO4--";
    bip_value = 9.4905;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "CO3--";
    bip_value = 3.3516;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "HCO3-";
    bip_value = -0.38795;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    // CO2

    species_name_1 = "CO2(aq)";
    set_bip(uij_T, id_map, species_name_1, species_name_1, 0.3587);

    species_name_2 = "H+";
    bip_value = 0.0;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Na+";
    bip_value = -0.436;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Ca++";
    bip_value = 0.0;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Ba++";
    bip_value = 0.0;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Sr++";
    bip_value = 0.0;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Mg++";
    bip_value = -2.855;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "OH-";
    bip_value = 0.0;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Cl-";
    bip_value = 15.015;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "SO4--";
    bip_value = 4.7896;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "CO3--";
    bip_value = 0.0;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "HCO3-";
    bip_value = -3.734;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    // H+

    species_name_1 = "H+";
    set_bip(uij_T, id_map, species_name_1, species_name_1, 0.0);

    species_name_2 = "Na+";
    bip_value = 0.0;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Ca++";
    bip_value = 0.0;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Ba++";
    bip_value = 0.0;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Sr++";
    bip_value = 0.0;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Mg++";
    bip_value = 0.0;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "OH-";
    bip_value = 0.0;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Cl-";
    bip_value = 0.0;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "SO4--";
    bip_value = 0.0;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "CO3--";
    bip_value = 0.0;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "HCO3-";
    bip_value = 0.0;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    // Na+

    species_name_1 = "Na+";
    set_bip(uij_T, id_map, species_name_1, species_name_1, 0.0);

    species_name_2 = "Ca++";
    bip_value = -4.656;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Ba++";
    bip_value = 2.338;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Sr++";
    bip_value = -0.62;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Mg++";
    bip_value = 1.3394;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "OH-";
    bip_value = 20.278;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Cl-";
    bip_value = 15.635;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "SO4--";
    bip_value = 11.681;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "CO3--";
    bip_value = 3.782;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "HCO3-";
    bip_value = 1.829;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    // Ca++

    species_name_1 = "Ca++";
    set_bip(uij_T, id_map, species_name_1, species_name_1, 0.0);

    species_name_2 = "Ba++";
    bip_value = 0.0;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Sr++";
    bip_value = -4.2533;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Mg++";
    bip_value = 5.1921;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "OH-";
    bip_value = 3.6084;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Cl-";
    bip_value = 11.14;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "SO4--";
    bip_value = 50.446;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "CO3--";
    bip_value = 0.0;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "HCO3-";
    bip_value = 0.0;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    // Ba++

    species_name_1 = "Ba++";
    set_bip(uij_T, id_map, species_name_1, species_name_1, 0.0);

    species_name_2 = "Sr++";
    bip_value = 0.0;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Mg++";
    bip_value = 0.0;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "OH-";
    bip_value = 0.0;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Cl-";
    bip_value = 14.89;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "SO4--";
    bip_value = 0.0;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "CO3--";
    bip_value = 0.0;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "HCO3-";
    bip_value = 0.0;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    // Sr++

    species_name_1 = "Sr++";
    set_bip(uij_T, id_map, species_name_1, species_name_1, 0.0);

    species_name_2 = "Mg++";
    bip_value = -1.437;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "OH-";
    bip_value = 0.0;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Cl-";
    bip_value = 15.689;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "SO4--";
    bip_value = 10;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "CO3--";
    bip_value = 0.0;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "HCO3-";
    bip_value = 0.0;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    // Mg++

    species_name_1 = "Mg++";
    set_bip(uij_T, id_map, species_name_1, species_name_1, 0.0);

    species_name_2 = "OH-";
    bip_value = 0.0;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "Cl-";
    bip_value = 12.132;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "SO4--";
    bip_value = 2.2791;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "CO3--";
    bip_value = 1.0;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "HCO3-";
    bip_value = 1,0;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    // OH-

    species_name_1 = "OH-";
    set_bip(uij_T, id_map, species_name_1, species_name_1, 5.6169);

    species_name_2 = "Cl-";
    bip_value = 13.628;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "SO4--";
    bip_value = 8.5902;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "CO3--";
    bip_value = 2.7496;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "HCO3-";
    bip_value = 0.0;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    // Cl-

    species_name_1 = "Cl-";
    set_bip(uij_T, id_map, species_name_1, species_name_1, 14.436);

    species_name_2 = "SO4--";
    bip_value = 12.407;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "CO3--";
    bip_value = 5.7267;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "HCO3-";
    bip_value = 14.035;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    // SO4--

    species_name_1 = "SO4--";
    set_bip(uij_T, id_map, species_name_1, species_name_1, 8.3194);

    species_name_2 = "CO3--";
    bip_value = 7.0067;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    species_name_2 = "HCO3-";
    bip_value = 6.9646;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    // CO3--

    species_name_1 = "CO3--";
    set_bip(uij_T, id_map, species_name_1, species_name_1, -1.3448);

    species_name_2 = "HCO3-";
    bip_value = 1.7241;
    set_bip(uij_T, id_map, species_name_1, species_name_2, bip_value);

    // HCO3-

    species_name_1 = "HCO3-";
    set_bip(uij_T, id_map, species_name_1, species_name_1, -0.0198);

    pimpl->linear_coeff_bips = uij_T;
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

auto EUNIQUACParams::ri(const std::string& name, double value) -> void
{
    pimpl->ri_values[name] = value;
}

auto EUNIQUACParams::ri(const std::map<std::string, double>& pairs) -> void
{
    for(const auto& pair : pairs)
        ri(pair.first, pair.second);
}

auto EUNIQUACParams::ri() const -> std::map<std::string, double>
{
    return pimpl->ri_values;
}

auto EUNIQUACParams::qi(const std::string& name, double value) -> void
{
    pimpl->qi_values[name] = value;
}

auto EUNIQUACParams::qi(const std::map<std::string, double>& pairs) -> void
{
    for(const auto& pair : pairs)
        qi(pair.first, pair.second);
}

auto EUNIQUACParams::qi() const -> std::map<std::string, double>
{
    return pimpl->qi_values;
}

auto EUNIQUACParams::set_uij_bips(
    const MatrixXd& uij_0_values,
    const MatrixXd& uij_T_values,
    const std::map<std::string, int>& species_id_map) -> void
{
    sanityCheckEnergyInteractionParams(uij_0_values, uij_T_values, species_id_map);
    pimpl->bips_species_id_map = species_id_map;
    pimpl->constant_coeff_bips = uij_0_values;
    pimpl->linear_coeff_bips = uij_T_values;

    // Guarantee that only the provided species in uij BIPs are used in qi and ri parameters
    std::map<std::string, double> new_qi_map, new_ri_map;
    const auto& old_qi_map = pimpl->qi_values;
    const auto& old_ri_map = pimpl->ri_values;
    for (const auto& species_and_id : species_id_map)
    {
        const auto& species_name = species_and_id.first;
        Assert(old_qi_map.find(species_name) != old_qi_map.end(),
            "The species set for u_ij BIPs has unknown E-UNIQUAC q_i parameters.",
            "One of the species while setting u_ij BIPs is lacking the q_i value.");

        Assert(old_ri_map.find(species_name) != old_ri_map.end(),
            "The species set for u_ij BIPs has unknown E-UNIQUAC r_i parameters.",
            "One of the species while setting u_ij BIPs is lacking the r_i value.");

        auto old_qi_value = old_qi_map.at(species_name);
        new_qi_map[species_name] = old_qi_value;
        auto old_ri_value = old_ri_map.at(species_name);
        new_ri_map[species_name] = old_ri_value;
    }

    pimpl->qi_values = new_qi_map;
    pimpl->ri_values = new_ri_map;
}

auto EUNIQUACParams::uij_0() const -> MatrixXd
{
    return pimpl->constant_coeff_bips;
}

auto EUNIQUACParams::uij_T() const -> MatrixXd
{
    return pimpl->linear_coeff_bips;
}

auto EUNIQUACParams::bips_species_id_map() const -> std::map<std::string, int>
{
    return pimpl->bips_species_id_map;
}

auto EUNIQUACParams::uij_0(
    const std::string& first_species_name,
    const std::string& second_species_name,
    double value) -> void
{
    Assert(!pimpl->bips_species_id_map.empty(),
        "Energy BIPs uij_0 cannot be set.",
        "The species id map should be provided before the BIPs initialization.");

    // TODO: implement exception handling when species names are unavailable
    const auto first_species_id = pimpl->bips_species_id_map.at(first_species_name);
    const auto second_species_id = pimpl->bips_species_id_map.at(second_species_name);
    // Energy BIPs are symmetric in E-UNIQUAC
    pimpl->constant_coeff_bips(first_species_id, second_species_id) = value;
    pimpl->constant_coeff_bips(second_species_id, first_species_id) = value;
}

auto EUNIQUACParams::uij_T(
    const std::string& first_species_name,
    const std::string& second_species_name,
    double value) -> void
{
    Assert(!pimpl->bips_species_id_map.empty(),
        "Energy BIPs uij_T cannot be set.",
        "The species id map should be provided before the BIPs initialization.");

    // TODO: implement exception handling when species names are unavailable
    const auto first_species_id = pimpl->bips_species_id_map.at(first_species_name);
    const auto second_species_id = pimpl->bips_species_id_map.at(second_species_name);
    // Energy BIPs are symmetric in E-UNIQUAC
    pimpl->linear_coeff_bips(first_species_id, second_species_id) = value;
    pimpl->linear_coeff_bips(second_species_id, first_species_id) = value;
}

auto EUNIQUACParams::bips_species_id_map(const std::map<std::string, int>& species_id_map) -> void
{
    pimpl->bips_species_id_map = species_id_map;
}

auto EUNIQUACParams::setDebyeHuckelGenericParameterA() -> void
{
    pimpl->useGeneralDebyeHuckelParameterA = true;
}

auto EUNIQUACParams::addNewSpeciesParameters(
    const std::string& species_name,
    double qi_value,
    double ri_value,
    const std::map<std::string, double>& u_0_values,
    const std::map<std::string, double>& u_T_values) -> void
{
    auto& species_id_map = pimpl->bips_species_id_map;
    auto& qi_map = pimpl->qi_values;
    auto& ri_map = pimpl->ri_values;

    // Assert that the new species is not defined in the species map, qi and ri maps
    Assert(species_id_map.find(species_name) == species_id_map.end(),
        "The species E-UNIQUAC parameters are already defined.",
        "The function addNewSpeciesParameters() is exclusive to add species which do not have any E-UNIQUAC "
        "parameters defined.");

    Assert(qi_map.find(species_name) == qi_map.end(),
        "The species E-UNIQUAC parameters are already defined.",
        "The function addNewSpeciesParameters() is exclusive to add species which do not have any E-UNIQUAC "
        "parameters defined.");

    Assert(ri_map.find(species_name) == ri_map.end(),
        "The species E-UNIQUAC parameters are already defined.",
        "The function addNewSpeciesParameters() is exclusive to add species which do not have any E-UNIQUAC "
        "parameters defined.");

    auto num_species = (int) species_id_map.size();
    auto id_new_species = num_species;
    species_id_map[species_name] = id_new_species;
    qi_map[species_name] = qi_value;
    ri_map[species_name] = ri_value;

    // Now we expand the BIPs matrices by one row and one column, both uninitialized at the end of each
    // dimension.
    auto& uij_0_matrix = pimpl->constant_coeff_bips;
    auto& uij_T_matrix = pimpl->linear_coeff_bips;
    uij_0_matrix.conservativeResize(num_species + 1, num_species + 1);
    uij_T_matrix.conservativeResize(num_species + 1, num_species + 1);

    // Adding new energy BIPs values, with j the index of existing species already stored
    double default_u_0_value = 2500.0;
    double default_u_T_value = 0.0;
    for (const auto& species_j_name_and_id: species_id_map)
    {
        const auto& species_j_name = species_j_name_and_id.first;
        auto species_j_id = species_j_name_and_id.second;

        // Add new uij_0 value, if provided. Otherwise, set it as the default value.
        bool has_uij_0_value = u_0_values.find(species_j_name) != u_0_values.end();
        if (has_uij_0_value) {
            auto new_uij_0 = u_0_values.at(species_j_name);
            uij_0_matrix(id_new_species, species_j_id) = new_uij_0;
            uij_0_matrix(species_j_id, id_new_species) = new_uij_0;
        }
        else {
            if (id_new_species == species_j_id)
                uij_0_matrix(id_new_species, species_j_id) = 0.0;
            else {
                uij_0_matrix(id_new_species, species_j_id) = default_u_0_value;
                uij_0_matrix(species_j_id, id_new_species) = default_u_0_value;
            }
        }

        // Add new uij_T value, if provided. Otherwise, set it as the default value.
        bool has_uij_T_value = u_T_values.find(species_j_name) != u_T_values.end();
        if (has_uij_T_value) {
            auto new_uij_T = u_T_values.at(species_j_name);
            uij_T_matrix(id_new_species, species_j_id) = new_uij_T;
            uij_T_matrix(species_j_id, id_new_species) = new_uij_T;
        }
        else {
            uij_T_matrix(id_new_species, species_j_id) = default_u_T_value;
            uij_T_matrix(species_j_id, id_new_species) = default_u_T_value;
        }
    }
}

auto EUNIQUACParams::useDebyeHuckelGenericParameterA() const -> bool
{
    return pimpl->useGeneralDebyeHuckelParameterA;
}
}  // namespace Reaktoro
