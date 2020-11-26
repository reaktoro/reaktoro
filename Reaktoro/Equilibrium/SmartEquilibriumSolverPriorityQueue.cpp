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

#include "SmartEquilibriumSolverPriorityQueue.hpp"

// C++ includes
#include <algorithm>
#include <tuple>

namespace Reaktoro {


SmartEquilibriumSolverPriorityQueue::SmartEquilibriumSolverPriorityQueue(const ChemicalSystem &system) : SmartEquilibriumSolverBase(Partition(system))
{
}

SmartEquilibriumSolverPriorityQueue::SmartEquilibriumSolverPriorityQueue(const Partition& partition) : SmartEquilibriumSolverBase(partition)
{
}

SmartEquilibriumSolverPriorityQueue::~SmartEquilibriumSolverPriorityQueue()
{
}

auto SmartEquilibriumSolverPriorityQueue::learn(ChemicalState& state, double T, double P, VectorConstRef be) -> void
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
    _properties = solver.properties();

    _result.timing.learn_chemical_properties = toc(CHEMICAL_PROPERTIES_STEP);

    //---------------------------------------------------------------------
    // SENSITIVITY MATRIX COMPUTATION STEP DURING THE LEARNING PROCESS
    //---------------------------------------------------------------------
    tic(SENSITIVITY_STEP);

    const auto& sensitivity = solver.sensitivity();

    _result.timing.learn_sensitivity_matrix = toc(SENSITIVITY_STEP);

    //---------------------------------------------------------------------
    // ERROR CONTROL MATRICES ASSEMBLING STEP DURING THE LEARNING PROCESS
    //---------------------------------------------------------------------
    tic(ERROR_CONTROL_MATRICES)

    // The amounts of the species at the calculated equilibrium state
    n = state.speciesAmounts();

    // The amounts of the equilibrium species amounts at the calculated equilibrium state
    ne = n(ies);

    const auto& x = _properties.moleFractions();
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

    // Fetch chemical potentials and their derivatives
    u = _properties.chemicalPotentials();

    // Auxiliary references to the derivatives dn/db, dn/dT, dn/dP, and du/dn
    const auto& dndT = sensitivity.dndT;
    const auto& dndP = sensitivity.dndP;
    const auto& dndb = sensitivity.dndb;

    // Compute the matrices du/dT, du/dP, du/db
    dudT = u.ddn * dndT + u.ddT; // du/dT = ∂u/∂n*∂n/∂T + ∂u/∂T
    dudP = u.ddn * dndP + u.ddP; // du/dP = ∂u/∂n*∂n/∂P + ∂u/∂P
    dudb = u.ddn * dndb;         // du/du = ∂u/∂n*∂n/∂b


    // The vector u(iprimary) with chemical potentials of primary species
    um.noalias() = u.val(imajor);

    // The matrix du(imajor)/dbe with derivatives of chemical potentials (imajor species only)
    const auto dumdbe = dudb(imajor, iee);
    const auto dumdT = dudT(imajor, 0);
    const auto dumdP = dudP(imajor, 0);

    // Compute matrix Mbe = 1/up * dup/db
    Mbe.noalias() = diag(inv(um)) * dumdbe;
    MT.noalias()  = diag(inv(um)) * dumdT;
    MP.noalias()  = diag(inv(um)) * dumdP;

    _result.timing.learn_error_control_matrices = toc(ERROR_CONTROL_MATRICES);

    //---------------------------------------------------------------------
    // STORAGE STEP DURING THE LEARNING PROCESS
    //---------------------------------------------------------------------
    tic(STORAGE_STEP)

    // Store the computed solution into the knowledge database
    database.push_back({T, P, be, state, _properties, solver.sensitivity(), Mbe, MT, MP, imajor});

    // Update the priority queue
    // ----------------------------------------------
    // Add new element in the priority queue
    priority.push_back(priority.size());

    // Set its rank to zero
    ranking.push_back(0);

    _result.timing.learn_storage = toc(STORAGE_STEP);
}

auto SmartEquilibriumSolverPriorityQueue::estimate(ChemicalState& state, double T, double P, VectorConstRef be) -> void
{
    _result.estimate.accepted = false;

    // Skip estimation if no previous full computation has been done
    if(database.empty())
        return;

    // Relative and absolute tolerance parameters
    const auto reltol = options.reltol;

    //---------------------------------------------------------------------
    // SEARCH STEP DURING THE LEARNING PROCESS
    //---------------------------------------------------------------------
    tic(SEARCH_STEP)

    // Variations of the elements
    Vector dbe;
    double dT, dP;
    
    // Check if an entry in the database pass the error test.
    // It returns (`success`, `error`, `ispecies`), where
    //   - `success` is true if error test succeeds, false otherwise.
    //   - `error` is the first error violating the tolerance
    //   - `ispecies` is the index of the species that fails the error test
    auto pass_error_test = [&be, &dbe, &T, &dT, &P, &dP, &reltol](const Record& record) -> std::tuple<bool, double, Index>
    {
        using std::abs;
        using std::max;
        const auto& be0 = record.be;
        const auto& T0 = record.T;
        const auto& P0 = record.P;
        const auto& Mbe0 = record.Mb;
        const auto& MT0 = record.MT;
        const auto& MP0 = record.MP;
        const auto& imajor = record.imajor;

        dbe.noalias() = be - be0;
        dT = T - T0;
        dP = P - P0;

        double error = 0.0;
        const auto size = Mbe0.rows();
        for(auto i = 1; i <= size; ++i)
        {
            const double delta_mu = Mbe0.row(size - i).dot(dbe) + MT0[size - i] * dT + MP0[size - i] * dP;
            error = max(error, abs(delta_mu)); // start checking primary species with least amount first
            if(error >= reltol)
                return { false, error, size - i };
        }

        return { true, error, -1 };
    };

    auto irecord_prev = priority.begin();
    for(auto irecord=priority.begin(); irecord!=priority.end(); ++irecord)
    {
        const auto& record = database[*irecord];
        const auto& imajor = record.imajor;

        //---------------------------------------------------------------------
        // ERROR CONTROL STEP DURING THE LEARNING PROCESS
        //---------------------------------------------------------------------
        tic(ERROR_CONTROL_STEP)

        const auto [success, error, ispecies] = pass_error_test(record);

        if(success)
        {
            _result.timing.estimate_error_control = toc(ERROR_CONTROL_STEP);

            //---------------------------------------------------------------------
            // TAYLOR EXTRAPOLATION STEP DURING THE LEARNING PROCESS
            //---------------------------------------------------------------------
            tic(TAYLOR_STEP)

            const auto& be0 = record.be;
            const auto& P0 = record.state.pressure();
            const auto& T0 = record.state.temperature();
            const auto& n0 = record.state.speciesAmounts();
            const auto& dndb0 = record.sensitivity.dndb;
            const auto& dndP0 = record.sensitivity.dndP;
            const auto& dndT0 = record.sensitivity.dndT;

            // Fetch reference values restricted to equilibrium species only
            const auto& ne0 = n0(ies);
            const auto& dnedbe0 = dndb0(ies, iee);
            const auto& dnedbPe0 = dndP0(ies);
            const auto& dnedbTe0 = dndT0(ies);

            // Perform Taylor extrapolation
            ne.noalias() = ne0 + dnedbe0 * (be - be0) + dnedbPe0 * (P - P0) + dnedbTe0 * (T - T0);

            _result.timing.estimate_taylor = toc(TAYLOR_STEP);

            // Check if all projected species amounts are positive
            const double ne_min = min(ne);
            const double ne_sum = sum(ne);
            const auto eps_n = options.amount_fraction_cutoff * ne_sum;

            if(ne_min <= -eps_n)
                continue;

            _result.timing.estimate_search = toc(SEARCH_STEP);

            //---------------------------------------------------------------------
            // After the search is finished successfully
            //---------------------------------------------------------------------

            // Assign small values to all the amount  in the interval [cutoff, 0] (instead of mirroring above)
            for(unsigned int i = 0; i < ne.size(); ++i) if(ne[i] < 0) ne[i] = options.learning.epsilon;

            // Update the amounts of elements for the equilibrium species
            //state = node.state; // this line was removed because it was destroying kinetics simulations
            state.setSpeciesAmounts(ne, ies);

            // Update the chemical properties of the system
            _properties = record.properties;  // FIXME: We actually want to estimate props =properties0 + variation : THIS IS A TEMPORARY SOLUTION!!!

            //-----------------------------------------------------------------------
            // DATABASE PRIORITY AND RANKING UPDATE STEP DURING THE ESTIMATE PROCESS
            //-----------------------------------------------------------------------
            tic(PRIORITY_UPDATE_STEP)

            // Increase ranking of of the used record
            ranking[*irecord] += 1;

            // Make sure the indices in the priority are ordered such that:
            // rank[priority[i-1]] > rank[priority[i]]
            auto comp = [&](Index l, Index r) { return ranking[l] > ranking[r]; };
            if( !((irecord == priority.begin()) || (ranking[*irecord_prev] >= ranking[*irecord])) ) 
            {
                std::stable_sort(priority.begin(), irecord + 1, comp);
            }

            _result.timing.estimate_database_priority_update = toc(PRIORITY_UPDATE_STEP);

            //---------------------------------------------------------------------
            // After the search is finished successfully
            //---------------------------------------------------------------------

            // Assign small values to all the amount  in the interval [cutoff, 0] (instead of mirroring above)
            for(unsigned int i = 0; i < ne.size(); ++i) 
                if(ne[i] < 0) 
                    ne[i] = options.learning.epsilon;

            // Update the amounts of elements for the equilibrium species
            //state = record.state; // this line was removed because it was destroying kinetics simulations
            state.setSpeciesAmounts(ne, ies);

            // Make sure that pressure and temperature is set to the current one
            state.setTemperature(T);
            state.setPressure(P);

            // Update the chemical properties of the system
            _properties = record.properties;  // FIXME: We actually want to estimate props =properties0 + variation : THIS IS A TEMPORARY SOLUTION!!!
            
            _result.estimate.accepted = true;

            return;
        }
        else {
            irecord_prev = irecord;
            continue;
        }
    }

    _result.estimate.accepted = false;
}

auto SmartEquilibriumSolverPriorityQueue::outputInfo() const -> void
{
}

} // namespace Reaktoro

