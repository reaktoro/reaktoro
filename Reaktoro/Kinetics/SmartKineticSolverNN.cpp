// Reaktoro is a unified framework for modeling chemically reactive systems.
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

#include "SmartKineticSolverNN.hpp"

// C++ includes
#include <functional>

namespace Reaktoro {

/// Implementation of the SmartKineticSolverNN functionality

SmartKineticSolverNN::SmartKineticSolverNN(const ReactionSystem& reactions, const Partition& partition)
: SmartKineticSolverBase(reactions, partition)
{
}

SmartKineticSolverNN::~SmartKineticSolverNN()
{
}

auto SmartKineticSolverNN::learn(ChemicalState& state, double& t, double dt) -> void
{
    //---------------------------------------------------------------------
    // INITIALIZATION
    //---------------------------------------------------------------------

    // Initialize sensitivity matrix by the identity matrix
    benk_S.setIdentity();

    // Initialize the kinetics state with the data at times t0 and t0
    ODEState ode_state;
    ode_state.u0 = benk;
    ode_state.t0 = t;

    //---------------------------------------------------------------------
    // CONVENTIONAL TIME-INTEGRATION DURING THE LEARNING PROCESS
    //---------------------------------------------------------------------
    tic(INTEGRATE_STEP)

    // Compute `benk` by the conventional numerical integration
    ode.solve(t, dt, benk, benk_S);
    //ode.solve_implicit_1st_order(t, dt, benk, benk_S); // TODO: compare CVODE and 1st order implicit scheme

    _result.timing.learn_integration = toc(INTEGRATE_STEP);

    //---------------------------------------------------------------------
    // STORAGE STEP DURING THE LEARNING PROCESS
    //---------------------------------------------------------------------
    tic(STORAGE_STEP)

    // Save the result time (t0 + dt), the obtain elements and kinetic species amounts, and the sensitivity values
    ode_state.t = t;
    ode_state.u = benk;
    ode_state.dudu0 = benk_S;

    // Update the composition of the amounts of equilibrium elements and kinetic species
    be = benk.head(Ee);
    nk = benk.tail(Nk);
    state.setSpeciesAmounts(nk, iks);

    // Save the kinetic state to the tree of learned states
    database.emplace_back(KineticRecord{state, _properties, sensitivity, ode_state, rates});

    _result.timing.learn_storage = toc(STORAGE_STEP);

}

auto SmartKineticSolverNN::estimate(ChemicalState& state, double& t, double dt) -> void
{
// Skip estimation if no previous full computation has been done
    if(database.empty())
        return;

    // Define initial state of the problem
    Vector benk0 = benk;

    // Comparison function based on the Euclidean distance
    auto distancefn = [&](const KineticRecord& a, const KineticRecord& b)
    {
        // Fetch benk0 from each TreeNode saved as u0 in ChemicalState
        const auto& benk0_a = a.ode_state.u0;
        const auto& benk0_b = b.ode_state.u0;

        return (benk0_a - benk0).squaredNorm() < (benk0_b - benk0).squaredNorm();  // TODO: We need to extend this later with T and P contributions too
    };

    //---------------------------------------------------------------------
    // SEARCH STEP DURING THE ESTIMATE PROCESS
    //---------------------------------------------------------------------
    tic(SEARCH_STEP)

    // Find the reference element (nearest to the new state benk)
    auto it = std::min_element(database.begin(), database.end(), distancefn);

    _result.timing.estimate_search = toc(SEARCH_STEP);

    //---------------------------------------------------------------------
    // TAYLOR PREDICTION STEP DURING THE ESTIMATE PROCESS
    //---------------------------------------------------------------------
    tic(TAYLOR_STEP)

    // Fetch the data stored in the reference element
    const auto& benk0_ref = it->ode_state.u0;
    const auto& benk_ref = it->ode_state.u;
    const auto& dndn0_ref = it->ode_state.dudu0;

    // Algorithm:
    // the reference state contains:
    // u0 -> benk0_ref      is the initial condition of reference vector
    // u -> benk_ref        is already calculated by integration reference vector
    // du/du0 -> dndn0_ref  is the sensitivity w.r.t. the initial condition
    // new initial values  : benk0
    // the predicted state : benk_new = benk_ref + dndn0_ref * (benk0 - benk0_ref)

    // Perform smart estimation of benk
    Vector benk_new;
    benk_new.noalias() = benk_ref + dndn0_ref * (benk0 - benk0_ref);

    _result.timing.estimate_taylor = toc(TAYLOR_STEP);

    //---------------------------------------------------------------------
    // ERROR CONTROL STEP DURING THE ESTIMATE PROCESS
    //---------------------------------------------------------------------
    tic(ERROR_CONTROL_STEP)

    // Fetch the be and nk unknowns from vector benk = [be; nk]
    VectorConstRef be_new = benk_new.head(Ee);
    VectorConstRef nk_new = benk_new.tail(Nk);

    // Fetch properties of the reference state
    const auto& state_ref = it->chemical_state;
    const auto& properties_ref = it->properties;
    const auto& sensitivity_ref = it->sensitivity;
    const auto& rates_ref = it->rates;
    const auto& n_ref = state_ref.speciesAmounts();

    const auto& be_ref = benk_ref.head(Ee);
    const auto& nk_ref = benk_ref.tail(Nk);
    const auto& ne_ref = n_ref(ies);

    // Get the sensitivity derivatives dln(a) / dn
    MatrixConstRef dlnadn_ref = properties_ref.lnActivities().ddn; // TODO: this line is assuming all species are equilibrium specie! get the rows and columns corresponding to equilibrium species
    VectorConstRef lna_ref = properties_ref.lnActivities().val;
    MatrixConstRef dsdn_ref = sensitivity_ref.dndb;

    // Get the sensitivity derivatives w.r.t. the equilibrium species dln(a) / dne
    MatrixConstRef dlnaedne_ref = dlnadn_ref(ies, ies);
    VectorConstRef lnae_ref = lna_ref(ies);
    MatrixConstRef dsdne_ref = dsdn_ref(ies, iee);

    // Define vectors of equilibrium species ne, delta(ne), delta(lnae)
    Vector dne, dlnae;
    dne.noalias() = dsdne_ref * (be_new - be_ref); // delta(ne) = dn/db * (be - be0)
    ne.noalias() = ne_ref + dne;                   // ne = ne_ref + delta(ne)
    dlnae.noalias() = dlnaedne_ref * dne;          // delta(dlnae) = d (ln ae) / d(ne)

    // Fetch mole fractions
    const auto& x_ref = properties_ref.moleFractions().val;
    VectorConstRef xe_ref = x_ref(ies);
    VectorConstRef xk_ref = x_ref(iks);

    ///---------------------------------------------------------------------
    // ERROR CONTROL RELATED TO THE EQUILIBRIUM SPECIES
    //---------------------------------------------------------------------

    // Run the test checking the negative values of equilibrium species
    bool equilibrium_neg_amount_check = ne.minCoeff() > options.cutoff;
    if(!equilibrium_neg_amount_check)
        return;

    // Run the test checking the variations of ln(a) of equilibrium species
    auto pass_ln_activities_equilibrium_species_error_test = [&]() -> bool
    {
        // Loop through all the equilibrium species
        for(int i = 0; i < xe_ref.size(); ++i) {

            // If the fraction is too small, skip the variational check
            if (xe_ref[i] < options.mole_fraction_cutoff)
                continue;

            // Perform the variational check
            if (std::abs(dlnae[i]) > options.abstol + options.reltol * std::abs(lnae_ref[i])) {
                _result.estimate.failed_with_species = system.species(ies[i]).name();
                _result.estimate.failed_with_amount = ne[i];
                return false;

            }
        }
        return true;
    };
    bool equilibrium_variation_check = pass_ln_activities_equilibrium_species_error_test();
    if(!equilibrium_variation_check)
        return;

    // ---------------------------------------------------------------------
    // ERROR CONTROL RELATED TO THE KINETIC SPECIES
    // ---------------------------------------------------------------------

    // Run the test checking the variations rates of kinetics species
    auto pass_rate_variation_kinetic_species_error_test = [&]() -> bool
    {
        /// Auxiliary vectors
        Vector dnk;
        dnk.noalias() = nk_new - nk_ref;
        Vector dn;
        dn.resize(Nk + Ne);
        dn(ies) << dne;
        dn(iks) << dnk;

        // Initialize reaction rates
        Vector drates = rates_ref.ddn * dn;

        // Loop over all kinetics species
        for(Index i = 0; i < xk_ref.size(); ++i){
            // If the fraction is too small, skip the variational check
            if(xk_ref[i] < options.mole_fraction_cutoff)
                continue;
            if(std::abs(drates.array()[i]) > options.abstol + options.reltol * std::abs(rates_ref.val.array()[i]))
                return false;
        }
        return true;

    };
    bool kinetics_r_variation_check = pass_rate_variation_kinetic_species_error_test();
    if(!kinetics_r_variation_check)
        return;

    _result.timing.estimate_error_control = toc(ERROR_CONTROL_STEP);

    // Assign small values to all the amount in the interval [cutoff, 0] (instead of mirroring above)
    for(unsigned int i = 0; i < benk_new.size(); ++i)
        if(benk_new[i] < 0)
            benk_new[i] = options.equilibrium.epsilon;

    // Update the solution of kinetic problem by new estimated value
    benk = benk_new;

    // Update properties by the reference one
    _properties = properties_ref;  // FIXME: We actually want to estimate properties = properties0 + variation : THIS IS A TEMPORARY SOLUTION!!!
    _properties.update(T, P);

    // Mark estimated result as accepted
    _result.estimate.accepted = true;

    // Update the time
    t += dt;
}

auto SmartKineticSolverNN::outputInfo() const -> void
{
}

} // namespace Reaktoro
