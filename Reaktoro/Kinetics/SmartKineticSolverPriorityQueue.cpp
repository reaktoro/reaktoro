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

#include "SmartKineticSolverPriorityQueue.hpp"

// C++ includes
#include <functional>

namespace Reaktoro {

/// Implementation of the SmartKineticSolverPriorityQueue functionality

SmartKineticSolverPriorityQueue::SmartKineticSolverPriorityQueue(const ReactionSystem& reactions, const Partition& partition)
: SmartKineticSolverBase(reactions, partition)
{
}

SmartKineticSolverPriorityQueue::~SmartKineticSolverPriorityQueue()
{
}

auto SmartKineticSolverPriorityQueue::learn(ChemicalState& state, double& t, double dt) -> void
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
    // Calculate `benk` by the conventional numerical integration
    ode.solve(t, dt, benk, benk_S);
    //ode.solve_implicit_1st_order(t, dt, benk, benk_S); // TODO: compare CVODE and 1st order implicit scheme
    //ode.integrate(t, benk, t + dt, benk_S); // TODO: produces a delayed estimations, why?
    _result.timing.learn_integration = toc(INTEGRATE_STEP);

    // Save the sensitivity values, the result time, and the obtain species' amount
    ode_state.t = t;
    ode_state.u = benk;
    ode_state.dudu0 = benk_S;

    // Update the composition of the amounts of equilibrium elements and kinetic species
    be = benk.head(Ee);
    nk = benk.tail(Nk);
    state.setSpeciesAmounts(nk, iks);

    //---------------------------------------------------------------------
    // ERROR CONTROL MATRICES ASSEMBLING STEP DURING THE LEARNING PROCESS
    //---------------------------------------------------------------------
    tic(ERROR_CONTROL_MATRICES)

    // The amounts of the species at the calculated equilibrium state
    VectorConstRef n = state.speciesAmounts();

    // The amounts of the equilibrium species amounts at the calculated equilibrium state
    ne = n(ies);

    // Mole fraction of the species at the calculated equilibrium state
    const auto& x = _properties.moleFractions();

    // Mole fraction of the equilibrium species species at the calculated equilibrium state
    Vector xe = x.val(ies);

    // Update the canonical form of formula matrix Ae so that we can identify primary species
    canonicalizer.updateWithPriorityWeights(ne);

    // Set tolerances for species' elements and fractions
    const auto eps_n = options.amount_fraction_cutoff * sum(ne);
    const auto eps_x = options.mole_fraction_cutoff;

    // Get indices of major species 'imajor'
    auto imajorminor = canonicalizer.Q();
    const auto nummajor = std::partition(imajorminor.begin(), imajorminor.end(),
                                         [&](Index i) { return ne[i] >= eps_n && xe[i] >= eps_x; })
                          - imajorminor.begin();
    const auto imajor = imajorminor.head(nummajor);

    // The chemical potentials at the calculated equilibrium state
    u = _properties.chemicalPotentials();

    // Auxiliary references to the derivatives dn/db and du/dn
    const auto& dndb = sensitivity.dndb;
    const auto& dudn = u.ddn;

    // Compute the matrix du/db = du/dn * dn/db
    dudb = dudn * dndb;

    // The vector u(imajor) with chemical potentials of major species
    um.noalias() = u.val(imajor);

    // The vector u(imajor) with chemical potentials of major species
    um.noalias() = u.val(imajor);

    // The matrix du(imajor)/dbe with derivatives of chemical potentials (major species only)
    const auto dumdbe = dudb(imajor, iee);

    // Compute matrix Mbe = 1/um * dum/db
    Mbe.noalias() = diag(inv(um)) * dumdbe;

    _result.timing.learn_error_control_matrices = toc(ERROR_CONTROL_MATRICES);

    //---------------------------------------------------------------------
    // STORAGE STEP DURING THE LEARNING PROCESS
    //---------------------------------------------------------------------
    tic(STORAGE_STEP)

    // Save the kinetic state to the tree of learned states
    database.emplace_back(KineticRecord{state, _properties, sensitivity, Mbe, imajor, ode_state, rates});

    // Update the priority queue
    // ----------------------------------------------
    // Add new element in the priority queue
    kinetics_priority.push_back(kinetics_priority.size());

    // Set its rank to zero
    kinetics_ranking.push_back(0);

    _result.timing.learn_storage = toc(STORAGE_STEP);
}

auto SmartKineticSolverPriorityQueue::estimate(ChemicalState& state, double& t, double dt) -> void
{
    _result.estimate.accepted = false;

    // Skip estimation if no previous full computation has been done
    if(database.empty())
        return;

    // Define initial state of the problem
    Vector benk0 = benk;
    // Define initial state of the problem
    Vector be = benk.head(Ee);
    Vector nk = benk.tail(Nk);

    //---------------------------------------------------------------------
    // SEARCH STEP DURING THE ESTIMATE PROCESS
    //---------------------------------------------------------------------
    tic(SEARCH_STEP)

    // Variations of the elements
    Vector dbe;

    // Check if an entry in the database pass the error test.
    // It returns (`success`, `error`, `ispecies`), where
    //   - `success` is true if error test succeeds, false otherwise.
    //   - `error` is the first error violating the tolerance
    //   - `ispecies` is the index of the species that fails the error test
    auto pass_equilibrium_potential_error_test = [&](const auto& node) -> std::tuple<bool, double, Index>
    {
        using std::abs;
        using std::max;

        const auto& benk_ref = node.ode_state.u0;
        const auto& be_ref = benk_ref.head(Ee);
        const auto& Mb_ref = node.Mb;
        const auto& imajor_ref = node.imajor;

        dbe.noalias() = be - be_ref;

        double error = 0.0;
        for(auto i = 0; i < imajor_ref.size(); ++i) {
            error = max(error, abs(Mb_ref.row(i) * dbe));
            if(error >= options.tol)
                return {false, error, imajor_ref[i] };
        }

        return { true, error, -1 };
    };

    // Initialize the inode_prev node
    auto inode_prev = kinetics_priority.begin();

    // Loop through all the nodes in priority ques
    for(auto inode=kinetics_priority.begin(); inode != kinetics_priority.end(); ++inode)
    {
        const auto& node = database[*inode];

        //---------------------------------------------------------------------
        // SEARCH CONTROL DURING THE ESTIMATE PROCESS
        //---------------------------------------------------------------------
        tic(EQUILIBRIUM_SPECIES_ERROR_CONTROL_STEP)

        const auto [success, error, ispecies] = pass_equilibrium_potential_error_test(node);

        if(success)
        {
            _result.timing.estimate_error_control += toc(EQUILIBRIUM_SPECIES_ERROR_CONTROL_STEP);

            //---------------------------------------------------------------------
            // TAYLOR DURING THE ESTIMATE PROCESS
            //---------------------------------------------------------------------
            tic(TAYLOR_STEP)

            // Fetch the data stored in the reference element
            const auto& benk0_ref = node.ode_state.u0;
            const auto& benk_ref = node.ode_state.u;
            const auto& dndn0_ref = node.ode_state.dudu0;

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
            // ERROR CONTROL THE ESTIMATE PROCESS
            //---------------------------------------------------------------------
            tic(ERROR_CONTROL_STEP)

            // Fetch the be and nk unknowns from vector benk = [be; nk]
            VectorConstRef be_new = benk_new.head(Ee);
            VectorConstRef nk_new = benk_new.tail(Nk);

            // -------------------------------------------------------------------------------------------------------------
            // Check the variations of equilibrium species
            // -------------------------------------------------------------------------------------------------------------

            // Define the function checking the negativity of the equilibrium species amounts
            auto pass_negative_equilibrium_species_amounts_error_test = [&](const auto& node) -> std::tuple<bool, VectorConstRef>
            {

                // Fetch properties of the reference state
                const auto& state_ref = node.chemical_state;
                const auto& sensitivity_ref = node.sensitivity;

                const auto& N = ies.size() + iks.size();
                const auto& be_ref = benk_ref.head(Ee);
                const auto& n_ref = state_ref.speciesAmounts();
                const auto& ne_ref = n_ref(ies);

                // Define vectors of equilibrium species ne, delta(ne), delta(lnae)
                VectorConstRef dne = sensitivity_ref.dndb(ies, iee) * (be_new - be_ref); // delta(ne) = dn/db * (be - be0)
                ne.noalias() = ne_ref + dne;                                                                  // ne = ne_ref + delta(ne)

//                    // Check if all projected species amounts are positive
//                    const double ne_min = min(ne);
//                    const double ne_sum = sum(ne);
//                    const auto eps_n = options.amount_fraction_cutoff * ne_sum;
//                    return {ne_min > -eps_n, dne};

                return {ne.minCoeff() > options.cutoff, dne}; // this check results in smaller number of learnings

            };
            // Check if the `negative equilibrium species amounts` pass the test
            auto [is_neg_equilibrium_test_passed, dne] = pass_negative_equilibrium_species_amounts_error_test(node);
            if(!is_neg_equilibrium_test_passed)
                continue;

            // -------------------------------------------------------------------------------------------------------------
            // Check the variations of kinetic species
            // -------------------------------------------------------------------------------------------------------------

            // Define the function checking the variation of kinetics rates
            auto pass_kinetic_rate_variation_error_test = [&](const auto& node, VectorConstRef dne) -> bool
            {
                const auto& rates_ref = node.rates;
                const auto& nk_ref = benk_ref.tail(Nk);
                const auto& properties_ref = node.properties;

                // Initialize delta_n = [dne; dnk]
                Vector dnk;
                dnk.noalias() = nk_new - nk_ref;
                Vector dn;
                dn.resize(Nk + Ne);
                dn(ies) << dne;
                dn(iks) << dnk;

                // Initialize reaction rates
                Vector drates = rates_ref.ddn * dn;

                // Fetch mole fractions
                const auto& x_ref = properties_ref.moleFractions().val;
                VectorConstRef xk_ref = x_ref(iks);

                // Loop over kinetic species with significant fractions
                for(Index i = 0; i < xk_ref.size(); ++i){
                    // If the fraction is too small, skip the variational check
                    if(xk_ref[i] < options.mole_fraction_cutoff)
                        continue;
                    if(std::abs(drates.array()[i]) > options.abstol + options.reltol * std::abs(rates_ref.val.array()[i])){
                        return false;
                    }
                }
                return true;
            };
            // Check if the variation in the kinetics rates pass the test
            const auto is_kin_rate_variation_test_passed = pass_kinetic_rate_variation_error_test(node, dne);
            if(!is_kin_rate_variation_test_passed)
                continue;

            _result.timing.estimate_error_control += toc(ERROR_CONTROL_STEP);

            _result.timing.estimate_search = toc(SEARCH_STEP);

            // -------------------------------------------------------------------------------------------------------------
            // Update the ranking and priority
            // -------------------------------------------------------------------------------------------------------------

            // Increase the raking of the node
            kinetics_ranking[*inode] += 1;

            // Sort the priority queue after increasing the ranking of the node
            auto comp = [&](Index l, Index r) { return kinetics_ranking[l] > kinetics_ranking[r]; };
            if( !((inode == kinetics_priority.begin()) || (kinetics_ranking[*inode_prev] >= kinetics_ranking[*inode])) ) {
                std::stable_sort(kinetics_priority.begin(), inode + 1, comp);
            }

            // Assign small values to all the amount in the interval [cutoff, 0] (instead of mirroring above)
            for(unsigned int i = 0; i < benk_new.size(); ++i) if(benk_new[i] < 0) benk_new[i] = options.equilibrium.epsilon;

            // -------------------------------------------------------------------------------------------------------------
            // Update the solution of kinetic problem by new estimated value
            // -------------------------------------------------------------------------------------------------------------
            //state = node.chemical_state; // this update doesn't work good for kinetics
            benk = benk_new;

            // Update the chemical properties of the system
            //properties = node.properties;  // FIXME: We actually want to estimate properties = properties0 + variation : THIS IS A TEMPORARY SOLUTION!!!

            // Mark estimated result as accepted
            _result.estimate.accepted = true;

            // Update the time
            t += dt;

            return;
        }
        else {
            inode_prev = inode;
            continue;
        }
    }

    _result.estimate.accepted = false;

}

auto SmartKineticSolverPriorityQueue::outputInfo() const -> void
{
}

SmartKineticSolverPriorityQueuePrimary::SmartKineticSolverPriorityQueuePrimary(const ReactionSystem& reactions, const Partition& partition)
 : SmartKineticSolverBase(reactions, partition)
{
}

/// Implementation of the SmartKineticSolverPriorityQueuePrimary functionality

SmartKineticSolverPriorityQueuePrimary::~SmartKineticSolverPriorityQueuePrimary()
{
}

auto SmartKineticSolverPriorityQueuePrimary::learn(ChemicalState& state, double& t, double dt) -> void
{
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

    // Calculate `benk` by the conventional numerical integration
    ode.solve(t, dt, benk, benk_S);
    _result.timing.learn_integration = toc(INTEGRATE_STEP);

    // Save the sensitivity values, the result time, and the obtain species' amount
    ode_state.t = t;
    ode_state.u = benk;
    ode_state.dudu0 = benk_S;

    // Update the composition of the amounts of equilibrium elements and kinetic species
    be = benk.head(Ee);
    nk = benk.tail(Nk);
    state.setSpeciesAmounts(nk, iks);

    //---------------------------------------------------------------------
    // ERROR CONTROL MATRICES ASSEMBLING STEP DURING THE LEARNING PROCESS
    //---------------------------------------------------------------------
    tic(ERROR_CONTROL_MATRICES)

    // The amounts of the species at the calculated equilibrium state
    VectorConstRef n = state.speciesAmounts();

    // The amounts of the equilibrium species amounts at the calculated equilibrium state
    ne = n(ies);

    // Update the canonical form of formula matrix Ae so that we can identify primary species
    canonicalizer.updateWithPriorityWeights(ne);

    // The order of the equilibrium species as (primary, secondary)
    const auto& iorder = canonicalizer.Q();

    // Assemble the vector of indices of equilibrium species as (primary, secondary)
    VectorXi ips(ies.size());
    for(auto i = 0; i < ips.size(); ++i)  // TODO: type Indices should be alias to VectorXi, to avoid such kind of codes
        ips[i] = ies[iorder[i]];

    // The number of primary species among the equilibrium species (Np <= Ne)
    const auto& Np = canonicalizer.numBasicVariables();

    // Store the indices of primary and secondary species in state
    state.equilibrium().setIndicesEquilibriumSpecies(ips, Np);

    // The indices of the primary species at the calculated equilibrium state
    VectorXiConstRef iprimary = ips.head(Np);

    // The chemical potentials at the calculated equilibrium state
    u = _properties.chemicalPotentials();

    // Auxiliary references to the derivatives dn/db and du/dn
    const auto& dndb = equilibrium.sensitivity().dndb;
    const auto& dudn = u.ddn;

    // Compute the matrix du/db = du/dn * dn/db
    dudb = dudn * dndb;

    // The vector u(iprimary) with chemical potentials of primary species
    up.noalias() = u.val(iprimary);

    // The matrix du(iprimary)/dbe with derivatives of chemical potentials (primary species only)
    const auto dupdbe = dudb(iprimary, iee);

    // Compute matrix Mbe = 1/up * dup/db
    Mbe.noalias() = diag(inv(up)) * dupdbe;

    _result.timing.learn_error_control_matrices = toc(ERROR_CONTROL_MATRICES);

    //---------------------------------------------------------------------
    // STORAGE STEP DURING THE LEARNING PROCESS
    //---------------------------------------------------------------------
    tic(STORAGE_STEP)

    // Save the kinetic state to the tree of learned states
    database.emplace_back(KineticRecord{state, _properties, sensitivity, Mbe, iprimary, ode_state, rates});

    // Update the priority queue
    // ----------------------------------------------
    // Add new element in the priority queue
    kinetics_priority.push_back(kinetics_priority.size());

    // Set its rank to zero
    kinetics_ranking.push_back(0);

    _result.timing.learn_storage = toc(STORAGE_STEP);
}

auto SmartKineticSolverPriorityQueuePrimary::estimate(ChemicalState& state, double& t, double dt) -> void
{
    _result.estimate.accepted = false;

    // Skip estimation if no previous full computation has been done
    if(database.empty())
        return;

    // Relative and absolute tolerance parameters
    const auto reltol = options.tol;

    // Define initial state of the problem
    Vector benk0 = benk;
    // Define initial state of the problem
    Vector be = benk.head(Ee);
    Vector nk = benk.tail(Nk);

    //---------------------------------------------------------------------
    // SEARCH STEP DURING THE ESTIMATE PROCESS
    //---------------------------------------------------------------------
    tic(SEARCH_STEP)

    Vector dbe;

    // The function that checks if a record in the database pass the error test.
    // It returns (`success`, `error`, `iprimaryspecies`), where
    // * `success` is true if error test succeeds, false otherwise.
    // * `error` is the first error value violating the tolerance
    // * `iprimaryspecies` is the index of the primary species that fails the error test
    auto pass_equilibrium_potential_error_test = [&](const KineticRecord& node) -> std::tuple<bool, double, Index>
    {
        using std::abs;
        using std::max;
        const auto& state_ref = node.chemical_state;
        const auto& benk_ref = node.ode_state.u;
        const auto& be_ref = benk_ref.head(Ee);
        const auto& Mbe0 = node.Mb;
        const auto& isue_ref = state_ref.equilibrium().indicesStrictlyUnstableElements();

        dbe.noalias() = be - be_ref;

        double error = 0.0;
        const auto size = Mbe0.rows();
        for(auto i = 1; i <= size; ++i) {
            error = max(error, abs(Mbe0.row(size - i) * dbe)); // start checking primary species with least amount first
            if(error >= reltol)
                return { false, error, size - i };
        }

        return { true, error, -1 };
    };

    auto inode_prev = kinetics_priority.begin();
    for(auto inode=kinetics_priority.begin(); inode != kinetics_priority.end(); ++inode)
    {
        const auto& node = database[*inode];

        //---------------------------------------------------------------------
        // SEARCH CONTROL DURING THE ESTIMATE PROCESS
        //---------------------------------------------------------------------
        tic(EQUILIBRIUM_SPECIES_ERROR_CONTROL_STEP)

        const auto [success, error, ispecies] = pass_equilibrium_potential_error_test(node);

        if(success)
        {
            _result.timing.estimate_error_control += toc(EQUILIBRIUM_SPECIES_ERROR_CONTROL_STEP);

            //---------------------------------------------------------------------
            // TAYLOR DURING THE ESTIMATE PROCESS
            //---------------------------------------------------------------------
            tic(TAYLOR_STEP)

            // Fetch the data stored in the reference element
            const auto& benk0_ref = node.ode_state.u0;
            const auto& benk_ref = node.ode_state.u;
            const auto& dndn0_ref = node.ode_state.dudu0;

            // Algorithm:
            // the reference state : u, u0, S = du/du0, t, f = du/dt
            // new initial values  : u0_new
            // the predicted state : u_new = u + S * (u0_new - u0)

            // Clarification:
            // benk0 is new initial condition (u0_tilde)
            // it.state.u0    is the initial condition of reference vector (u0)
            // it.state.u    is already calculated by integration vector  (u)
            // it.state.dudu0 is the sensitivity w.r.t. the initial condition (S)

            // Perform smart estimation of benk
            Vector benk_new;
            benk_new.noalias() = benk_ref + dndn0_ref * (benk0 - benk0_ref);

            _result.timing.estimate_taylor = toc(TAYLOR_STEP);

            //---------------------------------------------------------------------
            // ERROR CONTROL THE ESTIMATE PROCESS
            //---------------------------------------------------------------------
            tic(ERROR_CONTROL_STEP)

            // Fetch the be and nk unknowns from vector benk = [be; nk]
            VectorConstRef be_new = benk_new.head(Ee);
            VectorConstRef nk_new = benk_new.tail(Nk);

            // -------------------------------------------------------------------------------------------------------------
            // Check the variations of equilibrium species
            // -------------------------------------------------------------------------------------------------------------

            auto pass_negative_equilibrium_species_amounts_error_test = [&](const auto& node) -> std::tuple<bool, VectorConstRef>
            {

                // Fetch properties of the reference state
                const auto& state_ref = node.chemical_state;
                const auto& sensitivity_ref = node.sensitivity;

                const auto& N = ies.size() + iks.size();
                const auto& be_ref = benk_ref.head(Ee);
                const auto& n_ref = state_ref.speciesAmounts();
                const auto& ne_ref = n_ref(ies);

                // Define vectors of equilibrium species ne, delta(ne), delta(lnae)
                VectorConstRef dne = sensitivity_ref.dndb(ies, iee) * (be_new - be_ref); // delta(ne) = dn/db * (be - be0)
                ne.noalias() = ne_ref + dne;                                                                  // ne = ne_ref + delta(ne)


                // Check if all projected species amounts are positive
                const double ne_min = min(ne);
                const double ne_sum = sum(ne);
                const auto eps_n = options.amount_fraction_cutoff * ne_sum;

                return {ne_min > -eps_n, dne};

                // Negative cutoff check for the equilibrium
                // ---------------------------------------------------
                //return {ne.minCoeff() > 1e-2 * options.cutoff, dne};
            };

            auto [is_neg_equilibrium_test_passed, dne] = pass_negative_equilibrium_species_amounts_error_test(node);

            if(!is_neg_equilibrium_test_passed)
                continue;

            // -------------------------------------------------------------------------------------------------------------
            // Check the variations of kinetic species
            // -------------------------------------------------------------------------------------------------------------

            auto pass_kinetic_rate_variation_error_test = [&](const auto& node, VectorConstRef dne) -> bool
            {

                const auto& rates_ref = node.rates;
                const auto& nk_ref = benk_ref.tail(Nk);
                const auto& properties_ref = node.properties;

                // Initialize delta_n = [dne; delta_nk]
                Vector dnk;
                dnk.noalias() = nk_new - nk_ref;
                Vector dn;
                dn.resize(Nk + Ne);
                dn(ies) << dne;
                dn(iks) << dnk;

                // Initialize reaction rates
                Vector drates = rates_ref.ddn * dn;

                // Fetch mole fractions
                const auto& x_ref = properties_ref.moleFractions().val;
                VectorConstRef xk_ref = x_ref(iks);

                for(Index i = 0; i < xk_ref.size(); ++i){
                    // If the fraction is too small, skip the variational check
                    if(xk_ref[i] < options.mole_fraction_cutoff)
                        continue;
                    if(std::abs(drates.array()[i]) > options.abstol + options.reltol * std::abs(rates_ref.val.array()[i])){
                        return false;
                    }
                }
                return true;
            };

            const auto is_kin_rate_variation_test_passed = pass_kinetic_rate_variation_error_test(node, dne);

            if(!is_kin_rate_variation_test_passed)
                continue;

            _result.timing.estimate_error_control += toc(ERROR_CONTROL_STEP);

            _result.timing.estimate_search = toc(SEARCH_STEP);
            // -------------------------------------------------------------------------------------------------------------
            // Update the ranking and priority
            // -------------------------------------------------------------------------------------------------------------

            kinetics_ranking[*inode] += 1;

            auto comp = [&](Index l, Index r) { return kinetics_ranking[l] > kinetics_ranking[r]; };
            if( !((inode == kinetics_priority.begin()) || (kinetics_ranking[*inode_prev] >= kinetics_ranking[*inode])) ) {
                std::stable_sort(kinetics_priority.begin(), inode + 1, comp);
            }

            // Assign small values to all the amount in the interval [cutoff, 0] (instead of mirroring above)
            for(unsigned int i = 0; i < benk_new.size(); ++i) if(benk_new[i] < 0) benk_new[i] = options.equilibrium.epsilon;

            // -------------------------------------------------------------------------------------------------------------
            // Update the solution of kinetic problem by new estimated value
            // -------------------------------------------------------------------------------------------------------------
            //state = node.chemical_state;
            benk = benk_new;
            state = node.chemical_state;

            // Update the chemical properties of the system
            //properties = node.properties;  // FIXME: We actually want to estimate properties = properties0 + variation : THIS IS A TEMPORARY SOLUTION!!!

            // Mark estimated result as accepted
            _result.estimate.accepted = true;

            // Update the time
            t += dt;

            return;
        }
        else {
            inode_prev = inode;
            continue;
        }
    }

    _result.estimate.accepted = false;
}

auto SmartKineticSolverPriorityQueuePrimary::outputInfo() const -> void
{
}

} // namespace Reaktoro
