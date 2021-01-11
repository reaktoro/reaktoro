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

#include "SmartEquilibriumSolverNN.hpp"

// C++ includes
#include <algorithm>
#include <tuple>

namespace Reaktoro {

SmartEquilibriumSolverNN::SmartEquilibriumSolverNN(const ChemicalSystem &system) : SmartEquilibriumSolverBase(Partition(system))
{
}

SmartEquilibriumSolverNN::SmartEquilibriumSolverNN(const Partition& partition) : SmartEquilibriumSolverBase(partition)
{
}

SmartEquilibriumSolverNN::~SmartEquilibriumSolverNN()
{
}

auto SmartEquilibriumSolverNN::learn(ChemicalState& state, double T, double P, VectorConstRef be) -> void
{
    //---------------------------------------------------------------------
    // GIBBS ENERGY MINIMIZATION CALCULATION DURING THE LEARNING PROCESS
    //---------------------------------------------------------------------
    tic(EQUILIBRIUM_STEP)

    // Calculate the equilibrium state using conventional Gibbs energy minimization approach
    // and store the result of the Gibbs energy minimization calculation performed during learning
    _result.learning.gibbs_energy_minimization = solver.solve(state, T, P, be);

    // Check if the EquilibriumSolver calculation failed, if so, use cold-start
    if(!_result.learning.gibbs_energy_minimization.optimum.succeeded)
    {
        state.setSpeciesAmounts(0.0);
        _result.learning.gibbs_energy_minimization = solver.solve(state, T, P, be);
        if(!_result.learning.gibbs_energy_minimization.optimum.succeeded)
                return;
        
    }

    _result.timing.learn_gibbs_energy_minimization = toc(EQUILIBRIUM_STEP);

    //---------------------------------------------------------------------
    // CHEMICAL PROPERTIES UPDATE STEP DURING THE LEARNING PROCESS
    //---------------------------------------------------------------------
    tic(CHEMICAL_PROPERTIES_STEP)

    // Update the chemical properties of the system
    _properties =solver.properties();

    _result.timing.learn_chemical_properties = toc(CHEMICAL_PROPERTIES_STEP);

    //---------------------------------------------------------------------
    // STORAGE STEP DURING THE LEARNING PROCESS
    //---------------------------------------------------------------------
    tic(STORAGE_STEP)

    // Store the computed solution into the knowledge tree
    tree.push_back({be, state, _properties, solver.sensitivity()});

    _result.timing.learn_storage = toc(STORAGE_STEP);
}

auto SmartEquilibriumSolverNN::estimate(ChemicalState& state, double T, double P, VectorConstRef be) -> void
{
    // Skip estimation if no previous full computation has been done
    if(tree.empty())
        return;

    //---------------------------------------------------------------------
    // SEARCH STEP DURING THE ESTIMATE PROCESS
    //---------------------------------------------------------------------
    tic(SEARCH_STEP)

    // Comparison function based on the Euclidean distance
    auto distancefn = [&](const TreeNode& a, const TreeNode& b)
    {
        const auto& be_a = a.be;
        const auto& be_b = b.be;
        return (be_a - be).squaredNorm() < (be_b - be).squaredNorm();  // TODO: We need to extend this later with T and P contributions too (Allan, 26.06.2019)
    };

    // Find the entry with minimum "input" distance
    auto record = std::min_element(tree.begin(), tree.end(), distancefn);

    _result.timing.estimate_search = toc(SEARCH_STEP);

    //---------------------------------------------------------------------
    // TAYLOR STEP DURING THE ESTIMATE PROCESS
    //---------------------------------------------------------------------
    tic(TAYLOR_STEP)

    // Get all the data stored in the reference element
    const auto& be0 = record->be;
    const auto& state0 = record->state;
    const auto& properties0 = record->properties;
    const auto& sensitivity0 = record->sensitivity;
    const auto& n0 = state0.speciesAmounts();
    const auto& ne0 = n0(ies);

    // Get the sensitivity derivatives dln(a) / dn (assuming all species are equilibrium species)
    MatrixConstRef dlna_dn = properties0.lnActivities().ddn;
    VectorConstRef lna0 = properties0.lnActivities().val;
    MatrixConstRef dsdn0 = sensitivity0.dndb;

    // Get the sensitivity derivatives w.r.t. the equilibrium species dln(ae) / dne
    MatrixConstRef dlnae_dne = dlna_dn(ies, ies);
    VectorConstRef lnae0 = lna0(ies);
    MatrixConstRef dsdne0 = dsdn0(ies, iee);

    /// Auxiliary vectors delta(n) and delta(lna) in estimate function version v0
    Vector dne, dlnae;
    dne.noalias() = dsdne0 * (be - be0); // delta(ne) = dne/db * (b - b0)
    ne.noalias() = ne0 + dne;                                                   // n = n0 + delta(n)
    dlnae.noalias() = dlnae_dne * dne;                                          // delta(ln(a)) = d(lna)/dn * delta(n)

    _result.timing.estimate_taylor = toc(TAYLOR_STEP);

    //---------------------------------------------------------------------
    // ERROR CONTROL STEP DURING THE ESTIMATE PROCESS
    //---------------------------------------------------------------------
    tic(ERROR_CONTROL_STEP)

    // Fetch mole fractions
    const auto& x0 = properties0.moleFractions().val;
    Vector xe0 = x0(ies);

    // Perform the check for the negative amounts
    const bool amount_check = ne.minCoeff() > options.cutoff;
    if(!amount_check)
    {
        _result.timing.estimate_error_control = toc(ERROR_CONTROL_STEP);
        return;
    }

    // Check in the loop mole fractions, negative amount check, and variations of the ln(a)
    for(Index i = 0; i < ne.size(); ++i)
    {
        // If the fraction is too small, skip the variational check
        if(xe0[i] < options.mole_fraction_cutoff)
            continue;

        // Perform the variational check
        if(std::abs(dlnae[i]) > options.abstol + options.reltol * std::abs(lnae0[i])) 
        {
            _result.estimate.failed_with_species = system.species(ies[i]).name();
            _result.estimate.failed_with_amount = ne[i];
            _result.timing.estimate_error_control = toc(ERROR_CONTROL_STEP);
            return;
        }
    }

    _result.timing.estimate_error_control = toc(ERROR_CONTROL_STEP);

    //---------------------------------------------------------------------
    // After the search is finished successfully
    //---------------------------------------------------------------------

    // Assign small values to all the amount  in the interval [cutoff, 0] (instead of mirroring above)
    for(unsigned int i = 0; i < ne.size(); ++i) 
        if(ne[i] < 0) 
            ne[i] = options.learning.epsilon;

    // Update equilibrium species
    state.setSpeciesAmounts(ne, ies);

    // Make sure that pressure and temperature is set to the current one
    state.setTemperature(T);
    state.setPressure(P);

    // Set the estimate accepted status to true
    _result.estimate.accepted = true;

}

auto SmartEquilibriumSolverNN::outputInfo() const -> void
{
}

} // namespace Reaktoro

