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

#include "SmartKineticSolver.hpp"

// C++ includes
#include <functional>
#include <deque>

using namespace std::placeholders;

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/ChemicalVector.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Profiling.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumSolver.hpp>
#include <Reaktoro/Kinetics/SmartKineticOptions.hpp>
#include <Reaktoro/Kinetics/SmartKineticResult.hpp>
#include <Reaktoro/Optimization/Canonicalizer.hpp> // class needed for to fetch primary species
#include <Reaktoro/ODML/ClusterConnectivity.hpp> // class needed to maintain connectivity of the clusters
#include <Reaktoro/ODML/PriorityQueue.hpp> // class needed to maintain priority queues of clusters and records inside them

namespace Reaktoro {

namespace detail {

/// Return the hash number of a vector.
/// https://stackoverflow.com/questions/20511347/a-good-hash-function-for-a-vector
template<typename Vec>
auto hash(Vec& vec) -> std::size_t
{
    using T = decltype(vec[0]);
    const std::hash<T> hasher;
    std::size_t seed = vec.size();
    for(const auto& i : vec)
        seed ^= hasher(i) + 0x9e3779b9 + (seed << 6u) + (seed >> 2u);
    return seed;
}

/// Sort partially sorted array
void insertion_sort_from_assigned_start(VectorXi &array, Vector values, int start)
{
    // We obtain:
    // `array` is array of indices (of species) that must be correctly shuffled in correspondence to the `values`
    // `values` is array of values (amount of species) corresponding to `indices` (in chemical state)
    // `start` is the index, from which `values` is not sorted and `array` must be shuffled in correspondence
    //
    // Example of algorithm functionality:
    //
    // array: 35      16     9    28   1    4   18  24   31          34      33
    // values 35701.7 4947.9 59.8 57.9 14.3 1.2 0.8 0.03 8.91465e-05 1459.75 560.257
    // start: 9
    //
    // 1) array[9] = 34 must be inserted between 16 and 9 because 4947.9 > 1459.75 > 59.8
    //
    // array: 35      16     34       9    28   1    4   18  24   31          33
    // values 35701.7 4947.9 1459.75  59.8 57.9 14.3 1.2 0.8 0.03 8.91465e-05 560.257
    // start : 10
    //
    // 2) array[10] = 33 must be inserted between 34 and 9 because 1459.75 > 560.257 > 59.8

    Index i_value; // the key value that must be inserted
    double amount_i_value;

    // Run starting from the first "element to be inserted" indexed by `start` till the last
    for(int i = start, head = 0; i < array.size() && head < start; head++)
    {
        // Save `array[i]` that must be inserted later in the correct place
        i_value = array[i];
        amount_i_value  = values[i];

        // If we have reached the point, where elem value can be inserted
        // i.e., the value corresponding to `i`, `values[i]`, is bigger than `values[head]`
        if(values[i] > values[head])
        {
            // Move up smaller values
            for(auto k = i; k > head; k--)
            {
                array[k] = array[k - 1];
                values[k] = values[k - 1];
            }
            // Insert elem into its proper location
            array[head] = i_value;
            values[head] = amount_i_value;

            // Change the position of the next inserted element
            i++;
            // Refresh the head index (in case next value, values[i++] is bigger)
            head = 0;
        }
        else
        {
            // Otherwise, continue marching along 'array' increasing head
            head++;
            // If `head` reached `i`, we haven't found a place to insert
            if(head == i)
            {
                head = 0;   // reset head index
                i++;        // move on to the next element to insert
            }
        }
        // Otherwise, continue marching along 'array' increasing head
    }
    return;
}

} // namespace detail

/// A structure used to store the record of the knowledge database containing input, output, and derivatives data.
struct KineticRecord
{
    /// The temperature of the equilibrium state (in units of K).
    double T;

    /// The pressure of the equilibrium state (in units of Pa).
    double P;

    /// The amounts of elements in the equilibrium state (in units of mol).
    Vector be;

    /// The calculated equilibrium state at `T`, `P`, `be`.
    ChemicalState state;

    /// The chemical properties at the calculated equilibrium state.
    ChemicalProperties properties;

    /// The sensitivity derivatives at the calculated equilibrium state.
    EquilibriumSensitivity sensitivity;

    /// The matrix used to compute relative change of chemical potentials due to change in `be`.
    Matrix Mbe;

    /// Calculated kinetic state, containing initial and final value of [be, nk], initial and final time, and sensitivities
    ODEState ode_state;

    // Chemical kinetic rates
    ChemicalVector rates;
};

/// A structure used to store the node of tree for smart kinetic calculations.
struct SmartKineticNode{

    /// The calculated equilibrium state at `T`, `P`, `be`.
    ChemicalState chemical_state;

    /// The chemical properties at the calculated equilibrium state.
    ChemicalProperties properties;

    /// The sensitivity derivatives at the calculated equilibrium state.
    EquilibriumSensitivity sensitivity;

    /// The matrix used to compute relative change of chemical potentials due to change in `be`.
    Matrix Mb;

    // Indices of the major species
    VectorXi imajor;

    /// Calculated kinetic state, containing initial and final value of [be, nk], initial and final time, and sensitivities
    ODEState ode_state;

    // Chemical kinetic rates
    ChemicalVector rates;
};

struct SmartKineticSolver::Impl
{
    /// The kinetically-controlled chemical reactions
    ReactionSystem reactions;

    /// The chemical system instance
    ChemicalSystem system;

    /// The partition of the species in the chemical system
    Partition partition;

    /// The options of the kinetic solver
    SmartKineticOptions options;

    /// The result of the smart kinetic solver
    SmartKineticResult result;

    // The solver for solving the equilibrium equations using classical approach
    EquilibriumSolver equilibrium;

    /// The solver for solving the equilibrium equations using smart on-demand learning algorithm
    SmartEquilibriumSolver smart_equilibrium;

    /// The sensitivity of the equilibrium state
    EquilibriumSensitivity sensitivity;

    /// The chemical properties of the chemical system
    ChemicalProperties properties;

    /// The ODE solver instance
    ODESolver ode;

    /// The canonicalizer used to determine primary and secondary species
    Canonicalizer canonicalizer;

    /// The indices of the equilibrium and kinetic species
    Indices ies, iks;

    /// The indices of the elements in the equilibrium and kinetic partition
    Indices iee, ike;

    /// The number of equilibrium, kinetic, and all species
    Index Ne, Nk, N;

    /// The number of elements in the equilibrium partition
    Index Ee, E;

    /// The formula matrix of the equilibrium species
    Matrix Ae, Ak;

    /// The stoichiometric matrix w.r.t. the equilibrium and kinetic species
    Matrix Se, Sk;

    /// The coefficient matrix `A` of the kinetic rates
    Matrix A;

    /// The coefficient matrix `B` of the source rates
    Matrix B;

    /// The temperature of the chemical system (in units of K)
    double T;

    /// The pressure of the chemical system (in units of Pa)
    double P;

    /// The chemical potentials at the calculated equilibrium state
    ChemicalVector u;

    /// The amounts of the species in the chemical system
    Vector n;

    /// The molar composition of the equilibrium and kinetic species
    Vector ne, nk;

    /// The molar abundance of the elements in the equilibrium species
    Vector be;

    /// The scaled molar abundance of the elements in the equilibrium species
    Vector be_bar;

    /// The combined vector of elemental molar abundance and composition of kinetic species [be nk]
    Vector benk;

    /// The storage for matrix du/db = du/dn * dn/db
    Matrix dudb;

    /// The storage for vector u(iprimary)
    Vector up;

    /// The storage for matrix Mbe = inv(u(iprimary)) * du(iprimary)/db.
    Matrix Mbe;

    /// The vector with the values of the reaction rates
    ChemicalVector rates;

    /// The partial derivatives of the reaction rates `r` w.r.t. to `be`, `ne`, `nk`, `and `u = [be nk]`
    Matrix drdbe, drdne, drdnk, drdu;

    /// The source term
    ChemicalVector q;

    /// The partial derivatives of the source rates `q` w.r.t. to `be`, `ne`, `nk`, `and `u = [be nk]`
    Matrix dqdbe, dqdne, dqdnk, dqdu;

    /// The function that calculates the source term in the problem
    std::function<ChemicalVector(const ChemicalProperties&)> source_fn;

    /// The tree used to save the calculated equilibrium states and respective sensitivities
    std::deque<SmartKineticNode> tree;

    /// The sensitivity matrix of combined vector of elemental molar abundance and composition of kinetic species [be nk]
    Matrix benk_S;

    /// ------------------------------------------------
    /// Structures needed for the priority-based queue
    /// ------------------------------------------------

    /// The storage for vector u(imajor)
    Vector um;

    std::deque<Index> kinetics_priority;

    std::deque<Index> kinetics_ranking;

    /// ------------------------------------------------
    /// Structures needed for the priority-based cluster
    /// ------------------------------------------------

    /// The cluster storing learned input-output data with same classification.
    struct Cluster
    {
        /// The indices of the primary species for this cluster.
        VectorXi iprimary;

        /// The hash of the indices of the primary species for this cluster.
        std::size_t label = 0;

        /// The records stored in this cluster with learning data.
        std::deque<KineticRecord> records;

        /// The priority queue for the records based on their usage count.
        PriorityQueue priority;
    };

    /// The database containing the learned input-output data.
    struct Database
    {
        /// The clusters containing the learned input-output data points.
        std::deque<Cluster> clusters;

        /// The connectivity matrix of the clusters to determine how we move from one to another when searching.
        ClusterConnectivity connectivity;

        /// The priority queue for the clusters based on their usage counts.
        PriorityQueue priority;
    };

    /// The database with learned input-output data points.
    Database database;

    /// ---------------------------------------------------------------
    /// Structures needed for the priority-based cluster with extension
    /// ---------------------------------------------------------------

    /// The cluster (extended to kinetics) storing learned input-output data with same classification.
    struct ClusterExtended
    {
        /// The indices of the primary species for this cluster.
        VectorXi iprimary;

        /// The indices of the primary species and kinetics species for this cluster.
        VectorXi iprimary_ikin;

        /// The hash of the indices of the primary species for this cluster.
        std::size_t label = 0;

        /// The hash of the indices of the primary and kinetics species for this cluster.
        std::size_t label_extended = 0;

        /// The records stored in this cluster with learning data.
        std::deque<KineticRecord> records;

        /// The priority queue for the records based on their usage count.
        PriorityQueue priority;
    };

    /// The database containing the learned input-output data.
    struct DatabaseExtended
    {
        /// The clusters containing the learned input-output data points.
        std::deque<ClusterExtended> clusters;

        /// The connectivity matrix of the clusters to determine how we move from one to another when searching.
        ClusterConnectivity connectivity;

        /// The priority queue for the clusters based on their usage counts.
        PriorityQueue priority;
    };

    /// The database with learned input-output data points.
    DatabaseExtended database_extended;

    Impl(const ReactionSystem& reactions, const Partition& partition)
    : reactions(reactions),
      system(partition.system()),
      partition(partition),
      equilibrium(partition),
      smart_equilibrium(partition)
    {
        // Set the indices of the equilibrium, kinetic, and inert species
        ies = partition.indicesEquilibriumSpecies();
        iks = partition.indicesKineticSpecies();

        // Set the indices of the equilibrium and kinetic elements
        iee = partition.indicesEquilibriumElements();
        ike = partition.indicesKineticElements();

        // Set the number of equilibrium, kinetic, and inert species
        Ne = ies.size();
        Nk = iks.size();
        N = system.numSpecies();

        // Set the number of equilibrium and kinetic elements
        Ee = iee.size();
        E = system.numElements();

        // Initialise the formula matrix of the equilibrium partition
        Ae = partition.formulaMatrixEquilibriumPartition();

        // Initialize the canonicalizer with the formula matrix Ae of the equilibrium species
        canonicalizer.compute(Ae);

        // Initialise the formula matrix of the kinetic partition
        // Ak_tmp is of the size (Indices of elements that are included in kinetic species) x Nk
        Matrix Ak_tmp = partition.formulaMatrixKineticPartition();
        if(Nk)
        {
            // Initialize formula matrix of the kinetic partition
            // Ak is of the size E x Nk
            Ak = Matrix::Zero(E, Nk);

            for(Index i = 0; i < ike.size(); ++i)
            {
                // Copy the rows of Ak_tmp to the positions of the kinetic elements
                Ak.row(ike[i]) << Ak_tmp.row(i);
            }
        }

        // Initialise the stoichiometric matrices w.r.t. the equilibrium and kinetic species
        Se = cols(reactions.stoichiometricMatrix(), ies);
        Sk = cols(reactions.stoichiometricMatrix(), iks);

        // Initialise the coefficient matrix `A` of the kinetic rates
        A.resize(Ee + Nk, reactions.numReactions());
        A.topRows(Ee) = Ae * tr(Se);
        A.bottomRows(Nk) = tr(Sk);

        // Auxiliary identity matrix
        const Matrix I = identity(Ne + Nk, Ne + Nk);

        // Auxiliary selected equilibrium and kinetic rows of the identity matrix
        const Matrix Ie = rows(I, ies);
        const Matrix Ik = rows(I, iks);

        // Initialise the coefficient matrix `B` of the source rates
        B = zeros(Ee + Nk, N);
        B.topRows(Ee) = Ae * Ie;
        B.bottomRows(Nk) = Ik;

        // Allocate memory for the partial derivatives of the reaction rates `r` w.r.t. to `u = [be nk]`
        drdu.resize(reactions.numReactions(), Ee + Nk);

        // Allocate memory for the partial derivatives of the source rates `q` w.r.t. to `u = [be nk]`
        dqdu.resize(N, Ee + Nk);

        // Allocate the memory for the sensitivity matrix
        benk_S.resize(Ee + Nk, Ee + Nk);
    }

    auto setOptions(const SmartKineticOptions& _options) -> void
    {
        // Initialise the options of the kinetic solver
        this->options = _options;

        // Set options for the equilibrium calculations
        equilibrium.setOptions(options.equilibrium);
        smart_equilibrium.setOptions(options.smart_equilibrium);

    }

    auto addSource(ChemicalState state, double volumerate, const std::string& units) -> void
    {
        const Index num_species = system.numSpecies();
        const double volume = units::convert(volumerate, units, "m3/s");
        state.scaleVolume(volume);
        const Vector n = state.speciesAmounts();
        auto old_source_fn = source_fn;

        source_fn = [=](const ChemicalProperties& properties)
        {
            ChemicalVector q(num_species);
            q.val = n;
            if(old_source_fn)
                q += old_source_fn(properties);
            return q;
        };
    }

    auto addPhaseSink(const std::string& phase, double volumerate, const std::string& units) -> void
    {
        const double volume = units::convert(volumerate, units, "m3/s");
        const Index iphase = system.indexPhaseWithError(phase);
        const Index ifirst = system.indexFirstSpeciesInPhase(iphase);
        const Index size = system.numSpeciesInPhase(iphase);
        auto old_source_fn = source_fn;
        ChemicalScalar phasevolume;
        ChemicalVector q(size);

        source_fn = [=](const ChemicalProperties& properties) mutable
        {
            const auto n = properties.composition();
            const auto np = rows(n, ifirst, size);
            auto qp = rows(q, ifirst, size);
            phasevolume = properties.phaseVolumes()[iphase];
            qp = -volume*np/phasevolume;
            if(old_source_fn)
                q += old_source_fn(properties);
            return q;
        };
    }

    auto addFluidSink(double volumerate, const std::string& units) -> void
    {
        const double volume = units::convert(volumerate, units, "m3/s");
        const Indices& isolid_species = partition.indicesSolidSpecies();
        auto old_source_fn = source_fn;
        ChemicalScalar fluidvolume;
        ChemicalVector q;

        source_fn = [=](const ChemicalProperties& properties) mutable
        {
            const auto n = properties.composition();
            fluidvolume = properties.fluidVolume();
            q = -volume*n/fluidvolume;
            rows(q, isolid_species).fill(0.0);
            if(old_source_fn)
                q += old_source_fn(properties);
            return q;
        };
    }

    auto addSolidSink(double volumerate, const std::string& units) -> void
    {
        const double volume = units::convert(volumerate, units, "m3/s");
        const Indices& ifluid_species = partition.indicesFluidSpecies();
        auto old_source_fn = source_fn;
        ChemicalScalar solidvolume;
        ChemicalVector q;

        source_fn = [=](const ChemicalProperties& properties) mutable
        {
            const auto n = properties.composition();
            solidvolume = properties.solidVolume();
            q = -volume*n/solidvolume;
            rows(q, ifluid_species).fill(0.0);
            if(old_source_fn)
                q += old_source_fn(properties);
            return q;
        };
    }

    auto initialize(ChemicalState& state, double tstart) -> void {

        // Extract the composition of the equilibrium and kinetic species
        const auto &n = state.speciesAmounts();
        nk = n(iks);
        ne = n(ies);

        // Assemble the vector benk = [be nk]
        benk.resize(Ee + Nk);
        benk.head(Ee) = Ae * ne;
        benk.tail(Nk) = nk;

        // Initialize
        initialize(state, tstart, benk);

    }

    auto initialize(ChemicalState& state, double tstart, VectorConstRef benk) -> void
    {
        // Initialise the temperature and pressure variables
        T = state.temperature();
        P = state.pressure();

        // Define the ODE function
        ODEFunction ode_function = [&](double t, VectorConstRef u, VectorRef res)
        {
            return function(state, t, u, res);
        };

        // Define the jacobian of the ODE function
        ODEJacobian ode_jacobian = [&](double t, VectorConstRef u, MatrixRef res)
        {
            return jacobian(state, t, u, res);
        };

        // Initialise the ODE problem
        ODEProblem problem;
        problem.setNumEquations(Ee + Nk);
        problem.setFunction(ode_function);
        problem.setJacobian(ode_jacobian);

        // Define the options for the ODE solver
        ODEOptions options_ode = options.learning.ode;

        // Set the ODE problem and initialize the ODE solver
        ode.setProblem(problem);
        ode.setOptions(options_ode);
        ode.initialize(tstart, benk);

        // Set the options of the equilibrium solver
        smart_equilibrium.setOptions(options.learning.smart_equilibrium);
        equilibrium.setOptions(options.learning.equilibrium);

    }

    /// ------------------------------------------------------------------------------------------------------------ ///
    /// Learn functions
    /// ------------------------------------------------------------------------------------------------------------ ///

    /// Learn by performing a full kinetics calculation using integration of system of ODEs
    /// and save the learned state to the simple queue
    auto learn_nn_search(ChemicalState& state, double& t, double dt) -> void
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
        tic(INTEGRATE_STEP);
        // Compute `benk` by the conventional numerical integration
        ode.solve(t, dt, benk, benk_S);
        //ode.solve_implicit_1st_order(t, dt, benk, benk_S); // TODO: compare CVODE and 1st order implicit scheme
        //ode.integrate(t, benk, t + dt, benk_S); // TODO: produces a delayed estimations, why?
        result.timing.learn_integration = toc(INTEGRATE_STEP);

        //---------------------------------------------------------------------
        // STORAGE STEP DURING THE LEARNING PROCESS
        //---------------------------------------------------------------------
        tic(STORAGE_STEP);

        // Save the result time (t0 + dt), the obtain elements and kinetic species amounts, and the sensitivity values
        ode_state.t = t;
        ode_state.u = benk;
        ode_state.dudu0 = benk_S;

        // Update the composition of the amounts of equilibrium elements and kinetic species
        be = benk.head(Ee);
        nk = benk.tail(Nk);
        state.setSpeciesAmounts(nk, iks);

        // Save the kinetic state to the tree of learned states
        tree.emplace_back(SmartKineticNode{state, properties, sensitivity, Matrix::Zero(10, 10), VectorXi::Zero(10), ode_state, rates});

        result.timing.learn_storage = toc(STORAGE_STEP);

    }

    /// Learn by performing a full kinetics calculation using integration of system of ODEs
    /// and save the learned state (with the major potentials) to the priority queue
    auto learn_priority_based_acceptance_potential(ChemicalState& state, double& t, double dt) -> void
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
        tic(INTEGRATE_STEP);
        // Calculate `benk` by the conventional numerical integration
        ode.solve(t, dt, benk, benk_S);
        //ode.solve_implicit_1st_order(t, dt, benk, benk_S); // TODO: compare CVODE and 1st order implicit scheme
        //ode.integrate(t, benk, t + dt, benk_S); // TODO: produces a delayed estimations, why?
        result.timing.learn_integration = toc(INTEGRATE_STEP);

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
        tic(ERROR_CONTROL_MATRICES);

        // The amounts of the species at the calculated equilibrium state
        VectorConstRef n = state.speciesAmounts();

        // The amounts of the equilibrium species amounts at the calculated equilibrium state
        ne = n(ies);

        // Mole fraction of the species at the calculated equilibrium state
        const auto& x = properties.moleFractions();

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
        const auto u = properties.chemicalPotentials();

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

        result.timing.learn_error_control_matrices = toc(ERROR_CONTROL_MATRICES);

        //---------------------------------------------------------------------
        // STORAGE STEP DURING THE LEARNING PROCESS
        //---------------------------------------------------------------------
        tic(STORAGE_STEP);

        // Save the kinetic state to the tree of learned states
        tree.emplace_back(SmartKineticNode{state, properties, sensitivity, Mbe, imajor, ode_state, rates});

        // Update the priority queue
        // ----------------------------------------------
        // Add new element in the priority queue
        kinetics_priority.push_back(kinetics_priority.size());

        // Set its rank to zero
        kinetics_ranking.push_back(0);

        result.timing.learn_storage = toc(STORAGE_STEP);
    }

    /// Learn by performing a full kinetics calculation using integration of system of ODEs
    /// and save the learned state (with the primary potentials) to the priority queue
    auto learn_priority_based_acceptance_primary_potential(ChemicalState& state, double& t, double dt) -> void
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
        tic(INTEGRATE_STEP);
        // Calculate `benk` by the conventional numerical integration
        ode.solve(t, dt, benk, benk_S);
        result.timing.learn_integration = toc(INTEGRATE_STEP);

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
        tic(ERROR_CONTROL_MATRICES);

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
        const auto u = properties.chemicalPotentials();

        // Auxiliary references to the derivatives dn/db and du/dn
        const auto& dndb = equilibrium.sensitivity().dndb;
        const auto& dudn = u.ddn;

        // Compute the matrix du/db = du/dn * dn/db
        const Matrix dudb = dudn * dndb;

        // The vector u(iprimary) with chemical potentials of primary species
        const Vector up = u.val(iprimary);

        // The matrix du(iprimary)/dbe with derivatives of chemical potentials (primary species only)
        const auto dupdbe = dudb(iprimary, iee);

        // Compute matrix Mbe = 1/up * dup/db
        const Matrix Mbe = diag(inv(up)) * dupdbe;

        result.timing.learn_error_control_matrices = toc(ERROR_CONTROL_MATRICES);

        //---------------------------------------------------------------------
        // STORAGE STEP DURING THE LEARNING PROCESS
        //---------------------------------------------------------------------
        tic(STORAGE_STEP);

        // Save the kinetic state to the tree of learned states
        tree.emplace_back(SmartKineticNode{state, properties, sensitivity, Mbe, iprimary, ode_state, rates});

        // Update the priority queue
        // ----------------------------------------------
        // Add new element in the priority queue
        kinetics_priority.push_back(kinetics_priority.size());

        // Set its rank to zero
        kinetics_ranking.push_back(0);

        result.timing.learn_storage = toc(STORAGE_STEP);
    }

    /// Learn by performing a full kinetics calculation using integration of system of ODEs
    /// and save the learned state (with the primary potentials) to the priority-based cluster
    auto learn_clustering(ChemicalState& state, double& t, double dt) -> void
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
        tic(INTEGRATE_STEP);

        // Calculate `benk` by the conventional numerical integration
        ode.solve(t, dt, benk, benk_S);

        result.timing.learn_integration = toc(INTEGRATE_STEP);

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
        tic(ERROR_CONTROL_MATRICES);

        // The amounts of the species at the calculated equilibrium state
        Vector n = state.speciesAmounts();

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
        ChemicalVector u = properties.chemicalPotentials();

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

        result.timing.learn_error_control_matrices = toc(ERROR_CONTROL_MATRICES);

        //---------------------------------------------------------------------
        // STORAGE STEP DURING THE LEARNING PROCESS
        //---------------------------------------------------------------------
        tic(STORAGE_STEP);

        // Generate the hash number for the indices of primary species in the state
        const auto label = detail::hash(iprimary);

        // Find the index of the cluster that has same primary species
        auto iter = std::find_if(database.clusters.begin(), database.clusters.end(),
                                 [&](const Cluster& cluster) { return cluster.label == label; });

        // If cluster is found, store the new record in it, otherwise, create a new cluster
        if(iter < database.clusters.end())
        {
            auto& cluster = *iter;
            cluster.records.push_back({T, P, be, state, properties, sensitivity, Mbe, ode_state, rates});
            cluster.priority.extend();
        }
        else
        {
            // Create a new cluster
            Cluster cluster;
            cluster.iprimary = iprimary;
            cluster.label = label;
            cluster.records.push_back({T, P, be, state, properties, sensitivity, Mbe, ode_state, rates});
            cluster.priority.extend();

            // Append the new cluster in the database
            database.clusters.push_back(cluster);
            database.connectivity.extend();
            database.priority.extend();
        }

        result.timing.learn_storage = toc(STORAGE_STEP);
    }

    /// Learn by performing a full kinetics calculation using integration of system of ODEs
    /// and save the learned state (with the primary potentials) to the priority-based cluster
    auto learn_clustering_extended(ChemicalState& state, double& t, double dt) -> void
    {
        SmartKineticResult res = {};

        // Initialize sensitivity matrix by the identity matrix
        benk_S.setIdentity();

        // Initialize the kinetics state with the data at times t0 and t0
        ODEState ode_state;
        ode_state.u0 = benk;
        ode_state.t0 = t;

        //---------------------------------------------------------------------
        // CONVENTIONAL TIME-INTEGRATION DURING THE LEARNING PROCESS
        //---------------------------------------------------------------------
        tic(INTEGRATE_STEP);

        // Calculate `benk` by the conventional numerical integration
        ode.solve(t, dt, benk, benk_S);
        // Correct negative values by replacing them with small epsilon
        for(unsigned int i = 0; i < benk.size(); ++i) if(benk[i] < 0) benk[i] = options.equilibrium.epsilon;
        //ode.solve_implicit_1st_order(t, dt, benk, benk_S); // TODO: compare CVODE and 1st order implicit scheme
        //ode.integrate(t, benk, t + dt, benk_S); // TODO: produces a delayed estimations, why?
        result.timing.learn_integration += toc(INTEGRATE_STEP);

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
        tic(ERROR_CONTROL_MATRICES);

        // The amounts of the species at the calculated equilibrium state
        Vector n = state.speciesAmounts();

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

        // The indices of the primary species in the calculated equilibrium state
        VectorXiConstRef iprimary = ips.head(Np);

        //---------------------------------------------------------------------
        // CONSTRUCT AN EXTENDED IPRIMARY_IKIN VECTOR FOR CLUSTERING
        //---------------------------------------------------------------------

        // The indices of the primary and kinetic species the calculated equilibrium state
        VectorXi iprimary_ikin_tmp = VectorXi::Zero(iprimary.size() + iks.size());
        for (auto i = 0; i < iprimary.size(); ++i) iprimary_ikin_tmp[i] = iprimary[i];
        for (auto i = iprimary.size(); i < iprimary_ikin_tmp.size(); ++i) iprimary_ikin_tmp[i] = iks[i - iprimary.size()];

        // The amounts of species corresponding to indices iprimary_ikin
        Vector ne_iprimary_iks = Vector::Zero(iprimary_ikin_tmp.size());
        for (auto i = 0; i < iprimary_ikin_tmp.size(); ++i)
            ne_iprimary_iks[i] = state.speciesAmount(iprimary_ikin_tmp[i]);

        // Sort iprimary_ikin_tmp in correspondence with the values ne_iprimary_iks
        detail::insertion_sort_from_assigned_start(iprimary_ikin_tmp, ne_iprimary_iks, iprimary.size());

        // The indices of the primary and kinetic species the calculated equilibrium state (in correspondence with the values ne_iprimary_iks)
        VectorXiConstRef iprimary_ikin = iprimary_ikin_tmp;

        //---------------------------------------------------------------------

        // The chemical potentials at the calculated equilibrium state
        ChemicalVector u = properties.chemicalPotentials();

        // Auxiliary references to the derivatives dn/db and du/dn
        const auto& dndb = equilibrium.sensitivity().dndb;
        //const auto& dndb = options.use_smart_equilibrium_solver ? smart_equilibrium.sensitivity().dndb : equilibrium.sensitivity().dndb;
        const auto& dudn = u.ddn;

        // Compute the matrix du/db = du/dn * dn/db
        dudb = dudn * dndb;

        // The vector u(iprimary) with chemical potentials of primary species
        up.noalias() = u.val(iprimary);

        // The matrix du(iprimary)/dbe with derivatives of chemical potentials (primary species only)
        const auto dupdbe = dudb(iprimary, iee);

        // Compute matrix Mbe = 1/up * dup/db
        Mbe.noalias() = diag(inv(up)) * dupdbe;

        result.timing.learn_error_control_matrices = toc(ERROR_CONTROL_MATRICES);

        //---------------------------------------------------------------------
        // STORAGE STEP DURING THE LEARNING PROCESS
        //---------------------------------------------------------------------
        tic(STORAGE_STEP);

        // Generate the hash number for the indices of primary species in the state
        const auto label = detail::hash(iprimary);

        // Generate the hash number for the indices of primary and kinetic species in the state
        const auto label_extended = detail::hash(iprimary_ikin);

        // Find the index of the cluster that has same primary species
        auto iter = std::find_if(database_extended.clusters.begin(), database_extended.clusters.end(),
                                 [&](const ClusterExtended& cluster) { return cluster.label_extended == label_extended; });

        // If cluster is found, store the new record in it, otherwise, create a new cluster
        if(iter < database_extended.clusters.end())
        {
            auto& cluster = *iter;
            //cluster.records.push_back({T, P, be0, state, properties, sensitivity, Mbe, ode_state, rates});
            cluster.records.push_back({T, P, be, state, properties, sensitivity, Mbe, ode_state, rates});
            cluster.priority.extend();
        }
        else
        {
            // Create a new cluster
            ClusterExtended cluster;
            cluster.iprimary = iprimary;
            cluster.label = label;
            cluster.iprimary_ikin = iprimary_ikin;
            cluster.label_extended = label_extended;

            //cluster.records.push_back({T, P, be0, state, properties, sensitivity, Mbe, ode_state, rates});
            cluster.records.push_back({T, P, be, state, properties, sensitivity, Mbe, ode_state, rates});
            cluster.priority.extend();

            // Append the new cluster in the database
            database_extended.clusters.push_back(cluster);
            database_extended.connectivity.extend();
            database_extended.priority.extend();
        }

        result.timing.learn_storage = toc(STORAGE_STEP);
    }
    /// ------------------------------------------------------------------------------------------------------------ ///
    /// Estimate functions
    /// ------------------------------------------------------------------------------------------------------------ ///

    /// Estimate the equilibrium state using priority-based cluster for the search of the reference state and
    /// potentials of primary species for the acceptance criteria
    auto estimate_clustering_extended(ChemicalState& state, double& t, double dt) -> void
    {
        result.estimate.accepted = false;

        // Skip estimation if no previous full computation has been done
        if(database_extended.clusters.empty())
            return;

        // Relative and absolute tolerance parameters
        const auto reltol = options.tol;

        // Define initial state of the problem
        Vector benk0 = benk;
        // Define initial state of the problem
        Vector be = benk.head(Ee);
        Vector nk = benk.tail(Nk);

        const double T = state.temperature();
        const double P = state.pressure();

        // The threshold used to determine elements with insignificant amounts
        const auto eps_b = options.amount_fraction_cutoff * sum(be);

        // The current set of primary species in the chemical state
        const auto& iprimary = state.equilibrium().indicesPrimarySpecies();

        //---------------------------------------------------------------------
        // CONSTRUCT AN EXTENDED IPRIMARY_IKIN VECTOR FOR CLUSTERING
        //---------------------------------------------------------------------

        // The indices of the primary and kinetic species the calculated equilibrium state
        VectorXi iprimary_ikin_tmp = VectorXi::Zero(iprimary.size() + iks.size());
        for (auto i = 0; i < iprimary.size(); ++i) iprimary_ikin_tmp[i] = iprimary[i];
        for (auto i = iprimary.size(); i < iprimary_ikin_tmp.size(); ++i) iprimary_ikin_tmp[i] = iks[i - iprimary.size()];

        // The amounts of species corresponding to indices `iprimary_ikin`
        Vector ne_iprimary_iks = Vector::Zero(iprimary_ikin_tmp.size());
        for (auto i = 0; i < iprimary_ikin_tmp.size(); ++i)
            ne_iprimary_iks[i] = state.speciesAmount(iprimary_ikin_tmp[i]);

        // Sort `iprimary_ikin_tmp` in correspondence with the values `ne_iprimary_iks`
        // if `iprimary` is not empty, i.e., [35 16 9 28 18  4  1 24 31 | 34 33],
        // start sorting from the first element in the second part of the `iprimary_ikin_tmp` array
        // otherwise, [| 34 33], start from the second element
        if(iprimary.size())
            detail::insertion_sort_from_assigned_start(iprimary_ikin_tmp, ne_iprimary_iks, iprimary.size());
        else if (iks.size())
            detail::insertion_sort_from_assigned_start(iprimary_ikin_tmp, ne_iprimary_iks, 1);

        // The indices of the primary and kinetic species the calculated equilibrium state (in correspondence with the values ne_iprimary_iks)
        VectorXiConstRef iprimary_ikin = iprimary_ikin_tmp;
        //---------------------------------------------------------------------

        Vector dbe;

        // The function that checks if a record in the database pass the error test.
        // It returns (`success`, `error`, `iprimaryspecies`), where
        // * `success` is true if error test succeeds, false otherwise.
        // * `error` is the first error value violating the tolerance
        // * `iprimaryspecies` is the index of the primary species that fails the error test
        auto pass_error_test = [&be, &dbe, &reltol, &eps_b](const KineticRecord& record) -> std::tuple<bool, double, Index>
        {
            using std::abs;
            using std::max;
            const auto& state_ref = record.state;
            const auto& be_ref = record.be;
            const auto& Mbe_ref = record.Mbe;
            const auto& isue_ref = state_ref.equilibrium().indicesStrictlyUnstableElements();

            dbe.noalias() = be - be_ref;

            // Check if state_ref has strictly unstable elements (i.e. elements with zero amounts)
            // which cannot be used for Taylor estimation if positive amounts for those elements are given.
            if((dbe(isue_ref).array() > eps_b).any())
            {
                assert((be_ref(isue_ref).array() < eps_b).any()); // ensure this condition is never broken (effective during debug only)
                return { false, 9999, -1 };
            }

            double error = 0.0;
            const auto size = Mbe_ref.rows();
            for(auto i = 1; i <= size; ++i) {
                error = max(error, abs(Mbe_ref.row(size - i) * dbe)); // start checking primary species with least amount first
                if(error >= reltol)
                    return { false, error, size - i };
            }

            return { true, error, -1 };
        };

        // Generate the hash number for the indices of primary species in the state
        const auto label = detail::hash(iprimary);

        // Generate the hash number for the indices of primary species in the state
        const auto label_extended = detail::hash(iprimary_ikin);

        // The function that identifies the starting cluster index
        auto index_starting_cluster = [&]() -> Index
        {
            // If no primary species, then return number of clusters to trigger use of total usage counts of clusters
            if(iprimary.size() == 0)
                return database_extended.clusters.size(); // or database.clusters.size() - 1 ?

            // Find the index of the cluster with the same set of primary species (search those with highest count first)
            for(auto icluster : database_extended.priority.order())
                if(database_extended.clusters[icluster].label_extended == label_extended)
                    return icluster;

            // In no cluster with the same set of primary species if found, then return number of clusters
            return database_extended.clusters.size(); // or database.clusters.size() - 1 ?
        };

        // The index of the starting cluster
        const auto icluster = index_starting_cluster();

        // The ordering of the clusters to look for (starting with icluster)
        const auto& clusters_ordering = database_extended.connectivity.order(icluster);

        //---------------------------------------------------------------------
        // SEARCH STEP DURING THE ESTIMATE PROCESS
        //---------------------------------------------------------------------
        tic(SEARCH_STEP);

        int cluster_count = 0;
        // Iterate over all clusters (starting with icluster)
        for(auto jcluster : clusters_ordering)
        {
            // Fetch records from the cluster and the order they have to be processed in
            const auto& records = database_extended.clusters[jcluster].records;
            const auto& records_ordering = database_extended.clusters[jcluster].priority.order();

            int record_count = 0;
            // Iterate over all records in current cluster (using the  order based on the priorities)
            for(auto irecord : records_ordering)
            {
                const auto& record = records[irecord];

                //---------------------------------------------------------------------
                // ERROR CONTROL STEP DURING THE ESTIMATE PROCESS
                //---------------------------------------------------------------------
                tic(ERROR_CONTROL_STEP);

                // Check if the current record passes the error test
                const auto [success, error, iprimaryspecies] = pass_error_test(record);

                result.timing.estimate_error_control += toc(ERROR_CONTROL_STEP);

                if(success)
                {
                    //---------------------------------------------------------------------
                    // TAYLOR DURING THE ESTIMATE PROCESS
                    //---------------------------------------------------------------------
                    tic(TAYLOR_STEP);

                    // Fetch the data stored in the reference element
                    const auto& benk0_ref = record.ode_state.u0;
                    const auto& benk_ref = record.ode_state.u;
                    const auto& dndn0_ref = record.ode_state.dudu0;

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

                    result.timing.estimate_taylor = toc(TAYLOR_STEP);

                    //---------------------------------------------------------------------
                    // ERROR CONTROL THE ESTIMATE PROCESS
                    //---------------------------------------------------------------------
                    tic(ERROR_CONTROL_STEP);

                    // Fetch the be and nk unknowns from vector benk = [be; nk]
                    VectorConstRef be_new = benk_new.head(Ee);
                    VectorConstRef nk_new = benk_new.tail(Nk);

                    // -------------------------------------------------------------------------------------------------------------
                    // Check predicted negative values of equilibrium species
                    // -------------------------------------------------------------------------------------------------------------

                    // Define the function checking the negativity of the equilibrium species amounts
                    auto pass_negative_equilibrium_species_amounts_error_test = [&](const auto& record) -> std::tuple<bool, VectorConstRef, double, double>
                    {
                        // Fetch properties of the reference state
                        const auto& state_ref = record.state;
                        const auto& n_ref = state_ref.speciesAmounts();
                        const auto& P_ref = record.state.pressure();
                        const auto& T_ref = record.state.temperature();
                        const auto& dndb_ref = record.sensitivity.dndb;
                        const auto& dndP_ref = record.sensitivity.dndP;
                        const auto& dndT_ref = record.sensitivity.dndT;
                        const auto& be_ref = benk_ref.head(Ee);

                        const auto& ne_ref = n_ref(ies);
                        const auto& dnedbe_ref = dndb_ref(ies, iee);
                        const auto& dnedP_ref = dndP_ref(ies);
                        const auto& dnedT_ref = dndT_ref(ies);

                        // Define vectors of equilibrium species ne, delta(ne), delta(lnae)
                        VectorConstRef dne = dnedbe_ref * (be_new - be_ref) + dnedP_ref * (P - P_ref) + dnedT_ref * (T - T_ref);
                        ne.noalias() = ne_ref + dne;

                        // Check if all projected species amounts are positive
                        const double ne_min = min(ne);
                        const double ne_sum = sum(ne);
                        const auto eps_n = options.amount_fraction_cutoff * ne_sum;

                        return {ne_min > -eps_n, dne, ne_min, eps_n};
                        //return {ne.minCoeff() > options.cutoff, dne};
                    };
                    // Check if the `negative equilibrium species amounts` pass the test
                    auto [is_neg_equilibrium_test_passed, dne, ne_min, eps_n] = pass_negative_equilibrium_species_amounts_error_test(record);
                    if(!is_neg_equilibrium_test_passed){
                        continue;
                    }

                    // -------------------------------------------------------------------------------------------------------------
                    // Check the variations of kinetic species
                    // -------------------------------------------------------------------------------------------------------------

                    // Define the function checking the variations in the kinetics rate
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
                        Vector drates_kin = rates_ref.ddn(Eigen::all, iks) * dnk;

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
                    // Check if the `the variations in the kinetics rate` pass the test
                    const auto is_kin_rate_variation_test_passed = pass_kinetic_rate_variation_error_test(record, dne);
                    if(!is_kin_rate_variation_test_passed){
                        continue;
                    }

                    result.timing.estimate_error_control += toc(ERROR_CONTROL_STEP);

                    result.timing.estimate_search = toc(SEARCH_STEP);

                    // Assign small values to all the amount in the interval [cutoff, 0] (instead of mirroring above)
                    for(unsigned int i = 0; i < benk_new.size(); ++i) if(benk_new[i] < 0) benk_new[i] = options.equilibrium.epsilon;

                    // -------------------------------------------------------------------------------------------------------------
                    // Update the solution of kinetic problem by new estimated value
                    // -------------------------------------------------------------------------------------------------------------
                    benk = benk_new;
                    state = record.state; // ATTENTION: If this changes one day, make sure indices of equilibrium primary/secondary species, and indices of strictly unstable species/elements are also transfered from reference state to new state

                    // Update the chemical properties of the system
                    //properties = record.properties;  // FIXME: We actually want to estimate properties = properties0 + variation : THIS IS A TEMPORARY SOLUTION!!!

                    result.timing.estimate_taylor = toc(TAYLOR_STEP);

                    //---------------------------------------------------------------------
                    // DATABASE PRIORITY UPDATE STEP DURING THE ESTIMATE PROCESS
                    //---------------------------------------------------------------------
                    tic(PRIORITY_UPDATE_STEP);

                    // Increment priority of the current record (irecord) in the current cluster (jcluster)
                    database_extended.clusters[jcluster].priority.increment(irecord);

                    // Increment priority of the current cluster (jcluster) with respect to starting cluster (icluster)
                    database_extended.connectivity.increment(icluster, jcluster);

                    // Increment priority of the current cluster (jcluster)
                    database_extended.priority.increment(jcluster);

                    result.timing.estimate_database_priority_update = toc(PRIORITY_UPDATE_STEP);

                    // Mark estimated result as accepted
                    result.estimate.accepted = true;

                    // Update the time
                    t += dt;

                    return;
                }
                record_count++;
            }
            cluster_count++;
        }

        result.estimate.accepted = false;

    }


    /// Estimate the equilibrium state using priority-based cluster for the search of the reference state and
    /// potentials of primary species for the acceptance criteria
    auto estimate_clustering(ChemicalState& state, double& t, double dt) -> void
    {
        result.estimate.accepted = false;

        // Skip estimation if no previous full computation has been done
        if(database.clusters.empty())
            return;

        // Relative and absolute tolerance parameters
        const auto reltol = options.tol;

        // Define initial state of the problem
        Vector benk0 = benk;
        // Define initial state of the problem
        Vector be = benk.head(Ee);
        Vector nk = benk.tail(Nk);

        const double T = state.temperature();
        const double P = state.pressure();

        // The threshold used to determine elements with insignificant amounts
        const auto eps_b = options.amount_fraction_cutoff * sum(be);

        // The current set of primary species in the chemical state
        const auto& iprimary = state.equilibrium().indicesPrimarySpecies();

        Vector dbe;

        // The function that checks if a record in the database pass the error test.
        // It returns (`success`, `error`, `iprimaryspecies`), where
        // * `success` is true if error test succeeds, false otherwise.
        // * `error` is the first error value violating the tolerance
        // * `iprimaryspecies` is the index of the primary species that fails the error test
        auto pass_error_test = [&be, &dbe, &reltol, &eps_b](const KineticRecord& record) -> std::tuple<bool, double, Index>
        {
            using std::abs;
            using std::max;
            const auto& state_ref = record.state;
            const auto& be_ref = record.be;
            const auto& Mbe_ref = record.Mbe;
            const auto& isue_ref = state_ref.equilibrium().indicesStrictlyUnstableElements();

            dbe.noalias() = be - be_ref;

            // Check if state_ref has strictly unstable elements (i.e. elements with zero amounts)
            // which cannot be used for Taylor estimation if positive amounts for those elements are given.
            if((dbe(isue_ref).array() > eps_b).any())
            {
                assert((be_ref(isue_ref).array() < eps_b).any()); // ensure this condition is never broken (effective during debug only)
                return { false, 9999, -1 };
            }

            double error = 0.0;
            const auto size = Mbe_ref.rows();
            for(auto i = 1; i <= size; ++i) {
                error = max(error, abs(Mbe_ref.row(size - i) * dbe)); // start checking primary species with least amount first
                if(error >= reltol)
                    return { false, error, size - i };
            }

            return { true, error, -1 };
        };

        // Generate the hash number for the indices of primary species in the state
        const auto label = detail::hash(iprimary);

        // The function that identifies the starting cluster index
        auto index_starting_cluster = [&]() -> Index
        {
            // If no primary species, then return number of clusters to trigger use of total usage counts of clusters
            if(iprimary.size() == 0)
                return database.clusters.size(); // or database.clusters.size() - 1 ?

            // Find the index of the cluster with the same set of primary species (search those with highest count first)
            for(auto icluster : database.priority.order())
                if(database.clusters[icluster].label == label)
                    return icluster;

            // In no cluster with the same set of primary species if found, then return number of clusters
            return database.clusters.size(); // or database.clusters.size() - 1 ?
        };

        // The index of the starting cluster
        const auto icluster = index_starting_cluster();

        // The ordering of the clusters to look for (starting with icluster)
        const auto& clusters_ordering = database.connectivity.order(icluster);

        //---------------------------------------------------------------------
        // SEARCH STEP DURING THE ESTIMATE PROCESS
        //---------------------------------------------------------------------
        tic(SEARCH_STEP);

        // Iterate over all clusters (starting with icluster)
        for(auto jcluster : clusters_ordering)
        {
            // Fetch records from the cluster and the order they have to be processed in
            const auto& records = database.clusters[jcluster].records;
            const auto& records_ordering = database.clusters[jcluster].priority.order();

            // Iterate over all records in current cluster (using the  order based on the priorities)
            for(auto irecord : records_ordering)
            {
                const auto& record = records[irecord];

                //---------------------------------------------------------------------
                // ERROR CONTROL STEP DURING THE ESTIMATE PROCESS
                //---------------------------------------------------------------------
                tic(ERROR_CONTROL_STEP);

                // Check if the current record passes the error test
                const auto [success, error, iprimaryspecies] = pass_error_test(record);

                result.timing.estimate_error_control += toc(ERROR_CONTROL_STEP);

                if(success)
                {
                    //---------------------------------------------------------------------
                    // TAYLOR DURING THE ESTIMATE PROCESS
                    //---------------------------------------------------------------------
                    tic(TAYLOR_STEP);

                    // Fetch the data stored in the reference element
                    const auto& benk0_ref = record.ode_state.u0;
                    const auto& benk_ref = record.ode_state.u;
                    const auto& dndn0_ref = record.ode_state.dudu0;

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

                    result.timing.estimate_taylor = toc(TAYLOR_STEP);

                    //---------------------------------------------------------------------
                    // ERROR CONTROL THE ESTIMATE PROCESS
                    //---------------------------------------------------------------------
                    tic(ERROR_CONTROL_STEP);

                    // Fetch the be and nk unknowns from vector benk = [be; nk]
                    VectorConstRef be_new = benk_new.head(Ee);
                    VectorConstRef nk_new = benk_new.tail(Nk);

                    // -------------------------------------------------------------------------------------------------------------
                    // Check predicted negative values of equilibrium species
                    // -------------------------------------------------------------------------------------------------------------

                    // Define the function checking the negativity of the equilibrium species amounts
                    auto pass_negative_equilibrium_species_amounts_error_test = [&](const auto& record) -> std::tuple<bool, VectorConstRef>
                    {
                        // Fetch properties of the reference state
                        const auto& state_ref = record.state;
                        const auto& n_ref = state_ref.speciesAmounts();
                        const auto& P_ref = record.state.pressure();
                        const auto& T_ref = record.state.temperature();
                        const auto& dndb_ref = record.sensitivity.dndb;
                        const auto& dndP_ref = record.sensitivity.dndP;
                        const auto& dndT_ref = record.sensitivity.dndT;
                        const auto& be_ref = benk_ref.head(Ee);

                        const auto& ne_ref = n_ref(ies);
                        const auto& dnedbe_ref = dndb_ref(ies, iee);
                        const auto& dnedP_ref = dndP_ref(ies);
                        const auto& dnedT_ref = dndT_ref(ies);

                        // Define vectors of equilibrium species ne, delta(ne), delta(lnae)
                        VectorConstRef dne = dnedbe_ref * (be_new - be_ref) + dnedP_ref * (P - P_ref) + dnedT_ref * (T - T_ref);
                        ne.noalias() = ne_ref + dne;

                        // Check if all projected species amounts are positive
                        const double ne_min = min(ne);
                        const double ne_sum = sum(ne);
                        const auto eps_n = options.amount_fraction_cutoff * ne_sum;

                        return {ne_min > -eps_n, dne};
                    };
                    // Check if the `negative equilibrium species amounts` pass the test
                    auto [is_neg_equilibrium_test_passed, dne] = pass_negative_equilibrium_species_amounts_error_test(record);
                    if(!is_neg_equilibrium_test_passed)
                        continue;

                    // -------------------------------------------------------------------------------------------------------------
                    // Check the variations of kinetic species
                    // -------------------------------------------------------------------------------------------------------------

                    // Define the function checking the variations in the kinetics rate
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
                        //rates.val = rates_ref.val + drates;

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
                    // Check if the `the variations in the kinetics rate` pass the test
                    const auto is_kin_rate_variation_test_passed = pass_kinetic_rate_variation_error_test(record, dne);
                    if(!is_kin_rate_variation_test_passed)
                        continue;

                    result.timing.estimate_error_control += toc(ERROR_CONTROL_STEP);

                    result.timing.estimate_search = toc(SEARCH_STEP);

                    // Assign small values to all the amount in the interval [cutoff, 0] (instead of mirroring above)
                    for(unsigned int i = 0; i < benk_new.size(); ++i) if(benk_new[i] < 0) benk_new[i] = options.equilibrium.epsilon;

                    // -------------------------------------------------------------------------------------------------------------
                    // Update the solution of kinetic problem by new estimated value
                    // -------------------------------------------------------------------------------------------------------------
                    benk = benk_new;
                    state = record.state; // ATTENTION: If this changes one day, make sure indices of equilibrium primary/secondary species, and indices of strictly unstable species/elements are also transfered from reference state to new state

                    // Update the chemical properties of the system
                    //properties = record.properties;  // FIXME: We actually want to estimate properties = properties0 + variation : THIS IS A TEMPORARY SOLUTION!!!

                    result.timing.estimate_taylor = toc(TAYLOR_STEP);

                    //---------------------------------------------------------------------
                    // DATABASE PRIORITY UPDATE STEP DURING THE ESTIMATE PROCESS
                    //---------------------------------------------------------------------
                    tic(PRIORITY_UPDATE_STEP);

                    // Increment priority of the current record (irecord) in the current cluster (jcluster)
                    database.clusters[jcluster].priority.increment(irecord);

                    // Increment priority of the current cluster (jcluster) with respect to starting cluster (icluster)
                    database.connectivity.increment(icluster, jcluster);

                    // Increment priority of the current cluster (jcluster)
                    database.priority.increment(jcluster);

                    result.timing.estimate_database_priority_update = toc(PRIORITY_UPDATE_STEP);

                    // Mark estimated result as accepted
                    result.estimate.accepted = true;

                    // Update the time
                    t += dt;

                    return;
                }
            }
        }

        result.estimate.accepted = false;

    }

    /// Estimate the equilibrium state using the nearest neighbour for the search of the reference state and
    /// logarithms of activities for the acceptance criteria
    auto estimate_nn_search_acceptance_based_lna(ChemicalState& state, double& t, double dt) -> void
    {
        // Skip estimation if no previous full computation has been done
        if(tree.empty())
            return;

        // Define initial state of the problem
        Vector benk0 = benk;

        // Comparison function based on the Euclidean distance
        auto distancefn = [&](const SmartKineticNode& a, const SmartKineticNode& b)
        {
            // Fetch benk0 from each TreeNode saved as u0 in ChemicalState
            const auto& benk0_a = a.ode_state.u0;
            const auto& benk0_b = b.ode_state.u0;

            return (benk0_a - benk0).squaredNorm() < (benk0_b - benk0).squaredNorm();  // TODO: We need to extend this later with T and P contributions too
        };

        //---------------------------------------------------------------------
        // SEARCH STEP DURING THE ESTIMATE PROCESS
        //---------------------------------------------------------------------
        tic(SEARCH_STEP);

        // Find the reference element (nearest to the new state benk)
        auto it = std::min_element(tree.begin(), tree.end(), distancefn);

        result.timing.estimate_search = toc(SEARCH_STEP);

        //---------------------------------------------------------------------
        // TAYLOR PREDICTION STEP DURING THE ESTIMATE PROCESS
        //---------------------------------------------------------------------
        tic(TAYLOR_STEP);

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

        result.timing.estimate_taylor =  toc(TAYLOR_STEP);

        //---------------------------------------------------------------------
        // ERROR CONTROL STEP DURING THE ESTIMATE PROCESS
        //---------------------------------------------------------------------
        tic(ERROR_CONTROL_STEP);

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
                    result.estimate.failed_with_species = system.species(ies[i]).name();
                    result.estimate.failed_with_amount = ne[i];
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

        result.timing.estimate_error_control = toc(ERROR_CONTROL_STEP);

        // Assign small values to all the amount in the interval [cutoff, 0] (instead of mirroring above)
        for(unsigned int i = 0; i < benk_new.size(); ++i)
            if(benk_new[i] < 0)
                benk_new[i] = options.equilibrium.epsilon;

        // Update the solution of kinetic problem by new estimated value
        benk = benk_new;

        // Update properties by the reference one
        // properties = properties_ref; // does not influence results

        // Mark estimated result as accepted
        result.estimate.accepted = true;

        // Update the time
        t += dt;

    }

    /// Estimate the equilibrium state using the nearest neighbour for the search of the reference state and
    /// residual of the mass action equation for the acceptance criteria
    auto estimate_nn_search_acceptance_based_residual(ChemicalState& state, double& t, double dt) -> void
    {
        // Refresh acceptance result
        result.estimate.accepted = false;

        // Skip estimation if no previous full computation has been done
        if(tree.empty())
            return;

        // Define initial state of the problem
        Vector benk0 = benk;

        // Comparison function based on the Euclidean distance
        auto distancefn = [&](const SmartKineticNode& a, const SmartKineticNode& b)
        {
            // Fetch benk0 from each TreeNode saved as u0 in ChemicalState
            const auto& benk0_a = a.ode_state.u0;
            const auto& benk0_b = b.ode_state.u0;

            return (benk0_a - benk0).squaredNorm() < (benk0_b - benk0).squaredNorm();  // TODO: We need to extend this later with T and P contributions too
        };

        auto kinetics_reltol = options.reltol;
        auto kinetics_abstol = options.abstol;
        auto kinetics_cutoff = options.cutoff;
        auto fraction_tol = kinetics_abstol * 1e-2; // TODO: is it important to multiply by 1e-2

        //---------------------------------------------------------------------------------------
        // Step 1: Search for the reference element (closest to the new state input conditions)
        //---------------------------------------------------------------------------------------
        tic(SEARCH_STEP);

        // Find the reference element (nearest to the new state benk)
        auto it = std::min_element(tree.begin(), tree.end(), distancefn);

        result.timing.estimate_search = toc(SEARCH_STEP);

        //----------------------------------------------------------------------------
        // Step 2: Calculate predicted state with a first-order Taylor approximation
        //----------------------------------------------------------------------------
        tic(TAYLOR_STEP);

        // Fetch the data stored in the reference element
        const auto& benk0_ref = it->ode_state.u0;
        const auto& benk_ref = it->ode_state.u;
        const auto& dndn0_ref = it->ode_state.dudu0;

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

        result.timing.estimate_taylor = toc(TAYLOR_STEP);

        //----------------------------------------------
        // Step 3: Checking the acceptance criterion
        //----------------------------------------------
        tic(ERROR_CONTROL_STEP);

        // Fetch the be and nk unknowns from vector benk = [be; nk]
        VectorConstRef be_new = benk_new.head(Ee);
        VectorRef nk_new = benk_new.tail(Nk);

        // Fetch properties of the reference state
        const auto& state_ref = it->chemical_state;
        const auto& prop_ref = it->properties;
        const auto& sensitivity_ref = it->sensitivity;
        const auto& rates_ref = it->rates;

        // Get properties save in the chemical state
        const auto& T_ref = state_ref.temperature();
        const auto& P_ref = state_ref.pressure();
        const auto& y_ref = state_ref.equilibrium().elementChemicalPotentials();
        const auto& z_ref = state_ref.equilibrium().speciesStabilities();
        const auto u_ref = prop_ref.chemicalPotentials();
        const auto x_ref = prop_ref.moleFractions();

        // Get species amounts (equilibrium and kinetics) and element amounts
        const auto& N = ies.size() + iks.size();
        const auto& be_ref = benk_ref.head(Ee);
        const auto& nk_ref = benk_ref.tail(Nk);
        const auto& n_ref = state_ref.speciesAmounts();
        const auto& ne_ref = n_ref(ies);

        // Get species dual potentials, chemical potentials, and fractions related to the equilibrium species
        const auto& ze_ref = z_ref(ies);
        const auto& ue_ref = u_ref.val(ies);
        const auto& xe_ref = x_ref.val(ies);
        const auto& xk_ref = x_ref.val(iks);

        // Constant to scale the residual of the equilibrium species
        const auto RT = universalGasConstant*T_ref;

        /// Auxiliary vectors restricted to the equilibrium species
        Vector dne;
        Vector ze, ye, xe, ue, re;

        // Calculate perturbation of n
        dne.noalias() = sensitivity_ref.dndb(ies, iee) * (be_new - be_ref);

        ne.noalias() = ne_ref + dne;
        ze.noalias() = z_ref(ies);
        ye.noalias() = y_ref(iee);
        //std::cout << " y_ref =  " << y_ref.size() << ", y_ref        : " << tr(y_ref) << std::endl;

        /// -------------------------------------------------------------------------------------------------------------
        // Check if the predicted equilibrium species are not negative
        // -------------------------------------------------------------------------------------------------------------
        bool equilibrium_neg_amount_check = ne.minCoeff() > options.smart_equilibrium.cutoff;

//        // If cutoff test didn't pass, estimation has failded
//        if(equilibrium_neg_amount_check == false){
//            result.timing.estimate_error_control = toc(ERROR_CONTROL_STEP);
//            return;
//        }

        /// -------------------------------------------------------------------------------------------------------------
        // Check the variation of equilibrium species
        // -------------------------------------------------------------------------------------------------------------

        // Calculate u(bar) = u(ref) + dudT(ref)*dT + dudP(ref)*dP + dudn(ref)*dn
        ue = ue_ref + u_ref.ddn(ies, ies) *dne;

        // Calculate x(bar) = x(ref) + dxdn(ref)*dn + dxdP(ref)*dP + dxdn(ref)*dn
        xe = xe_ref + x_ref.ddn(ies, ies) * dne;

        // Obtain formula matrix of equilibrium partition
        MatrixConstRef Ae = partition.formulaMatrixEquilibriumPartition();

        //std::cout << " ye =  " << ye.size() << ", ye        : " << tr(ye) << std::endl;
        //getchar();
        // Calculate the equilibrium residuals of the equilibrium species
        re = abs(ue - tr(Ae)*ye - ze)/RT;  // TODO: We should actually collect the entries in u and z corresponding to equilibrium species
        //std::cout << " re =  " << re.size() << ", re        : " << tr(re) << std::endl;
        //getchar();

        // Eliminate species with mole fractions below cutoff from the residual analysis
        for(auto i = 0; i < ne.size(); ++i)
            if(xe[i] < options.smart_equilibrium.mole_fraction_cutoff)
                re[i] = 0.0; // set their residuals to zero

        // Estimate the residual error of the trial Taylor approximation for the equilibrium species
        Index ispecies;
        const double error = re.maxCoeff(&ispecies);
        bool equilibrium_variation_check = error <= options.tol;
        //std::cout << " re =  " << re.size() << ", re        : " << tr(re) << std::endl;
        //std::cout << " error =  " << error << std::endl;
        //std::cout << " options.tol =  " << options.tol << std::endl;
        //getchar();

//        // If equilibrium variation test didn't pass, estimation has failed
//        if(equilibrium_variation_check == false){
//            result.timing.estimate_error_control = toc(ERROR_CONTROL_STEP);
//            return;
//        }
        /// -------------------------------------------------------------------------------------------------------------
        // Check the variations of kinetic species
        // -------------------------------------------------------------------------------------------------------------

        Vector dnk;
        // Initialize delta_n = [delta_ne; delta_nk]
        dnk.noalias() = nk_new - nk_ref;

        Vector dn;
        dn.resize(N);
        dn(ies) << dne;
        dn(iks) << dnk;

        // Initialize reaction rates
        Vector drates = rates_ref.ddn * dn;
        //rates.val = rates_ref.val + drates;

        bool kinetics_r_variation_check = true;
        /// The vector of amounts of equilibrium species
        Vector xk;
        xk.noalias() = xk_ref + x_ref.ddn(iks, iks) * dnk;

        for(auto i = 0; i < nk.size(); ++i){
            if(xk[i] < options.mole_fraction_cutoff)
                continue; // if the fraction is too small, skip the variational check
            if(std::abs(drates.array()[i]) > options.abstol + options.reltol * std::abs(rates_ref.val.array()[i])){
                kinetics_r_variation_check = false;
                //result.timing.estimate_error_control = toc(ERROR_CONTROL_STEP);
                //return;
                break;
            }
        }
        result.timing.estimate_error_control = toc(ERROR_CONTROL_STEP);

        if(!equilibrium_variation_check || !equilibrium_neg_amount_check || !kinetics_r_variation_check)
            //if(!equilibrium_variation_check || !kinetics_r_variation_check)
        {
            /*
            if(step>=1000) {
                //std::cout << "|xk_ref|  : " << TR(xk_ref)  << std::endl;
                //std::cout << fraction_tol  : " << fraction_tol  << std::endl;
                //std::cout << "|delta_r| : " << delta_r.array().abs() << std::endl;
                //std::cout << "|r_ref|   : " << r_ref.val.array().abs() << std::endl;
                //std::cout << "abstol + reltol * |r_ref|    : " << kinetics_abstol + kinetics_reltol * std::abs(r_ref.val.array()[i]) << std::endl;

                std::cout << "lnae_var_check     : " << equilibrium_variation_check;
                std::cout << ", eq_amount_check  : " << equilibrium_neg_amount_check;
                std::cout << ", r_var_check    : " << kinetics_r_variation_check << std::endl;

                getchar();
            }
            */
            ///*
            //std::cout << "equilibrium_variation_check     : " << equilibrium_variation_check;
            //std::cout << ", equilibrium_neg_amount_check eq_amount_check  : " << equilibrium_neg_amount_check;
            //std::cout << ", kinetics_r_variation_check    : " << kinetics_r_variation_check << std::endl;
            //getchar();
            //*/
            return;
        }
        //----------------------------------------------
        // Step 3: Checking the acceptance criterion
        //----------------------------------------------

        // Update the chemical properties of the system
        //properties = prop_ref;  // FIXME: We actually want to estimate properties = properties0 + variation : THIS IS A TEMPORARY SOLUTION!!!

        for(unsigned int i = 0; i < benk_new.size(); ++i) if(benk_new[i] < 0) benk_new[i] = options.equilibrium.epsilon;

        // Update the solution of kinetic problem by new estimated value
        benk = benk_new;

        // Mark estimated result as accepted
        result.estimate.accepted = true;

        // Update the time
        t += dt;

    }

    /// Estimate the equilibrium state using the priority-based queue for the search of the reference state and
    /// potentials of major species for the acceptance criteria
    auto estimate_priority_based_acceptance_potential(ChemicalState& state, double& t, double dt) -> void
    {
        result.estimate.accepted = false;

        // Skip estimation if no previous full computation has been done
        if(tree.empty())
            return;

        // Relative and absolute tolerance parameters
        const auto reltol = options.tol;

        // Define initial state of the problem
        Vector benk0 = benk;
        // Define initial state of the problem
        Vector be = benk.head(Ee);;
        Vector nk = benk.tail(Nk);;

        //---------------------------------------------------------------------
        // SEARCH STEP DURING THE ESTIMATE PROCESS
        //---------------------------------------------------------------------
        tic(SEARCH_STEP);

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
            const auto& node = tree[*inode];

            //---------------------------------------------------------------------
            // SEARCH CONTROL DURING THE ESTIMATE PROCESS
            //---------------------------------------------------------------------
            tic(ERROR_CONTROL_STEP);

            const auto [success, error, ispecies] = pass_equilibrium_potential_error_test(node);

            if(success)
            {
                result.timing.estimate_error_control += toc(ERROR_CONTROL_STEP);

                //---------------------------------------------------------------------
                // TAYLOR DURING THE ESTIMATE PROCESS
                //---------------------------------------------------------------------
                tic(TAYLOR_STEP);

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

                result.timing.estimate_taylor = toc(TAYLOR_STEP);

                //---------------------------------------------------------------------
                // ERROR CONTROL THE ESTIMATE PROCESS
                //---------------------------------------------------------------------
                tic(ERROR_CONTROL_STEP);

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

                result.timing.estimate_error_control += toc(ERROR_CONTROL_STEP);

                result.timing.estimate_search = toc(SEARCH_STEP);

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
                result.estimate.accepted = true;

                // Update the time
                t += dt;

                return;
            }
            else {
                inode_prev = inode;
                continue;
            }
        }

        result.estimate.accepted = false;

    }

    /// Estimate the equilibrium state using the priority-based queue for the search of the reference state and
    /// potentials of primary species for the acceptance criteria
    auto estimate_priority_based_acceptance_primary_potential(ChemicalState& state, double& t, double dt) -> void
    {
        result.estimate.accepted = false;

        // Skip estimation if no previous full computation has been done
        if(tree.empty())
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
        tic(SEARCH_STEP);

        // The threshold used to determine elements with insignificant amounts
        const auto eps_b = options.amount_fraction_cutoff * sum(be);

        Vector dbe;

        // The function that checks if a record in the database pass the error test.
        // It returns (`success`, `error`, `iprimaryspecies`), where
        // * `success` is true if error test succeeds, false otherwise.
        // * `error` is the first error value violating the tolerance
        // * `iprimaryspecies` is the index of the primary species that fails the error test
        auto pass_equilibrium_potential_error_test = [&](const SmartKineticNode& node) -> std::tuple<bool, double, Index>
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
            const auto& node = tree[*inode];

            //---------------------------------------------------------------------
            // SEARCH CONTROL DURING THE ESTIMATE PROCESS
            //---------------------------------------------------------------------
            tic(ERROR_CONTROL_STEP);

            const auto [success, error, ispecies] = pass_equilibrium_potential_error_test(node);

            if(success)
            {
                result.timing.estimate_error_control += toc(ERROR_CONTROL_STEP);

                //---------------------------------------------------------------------
                // TAYLOR DURING THE ESTIMATE PROCESS
                //---------------------------------------------------------------------
                tic(TAYLOR_STEP);

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

                result.timing.estimate_taylor = toc(TAYLOR_STEP);

                //---------------------------------------------------------------------
                // ERROR CONTROL THE ESTIMATE PROCESS
                //---------------------------------------------------------------------
                tic(ERROR_CONTROL_STEP);

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
                //std::cout << "is_neg_equilibrium_test_passed = " << is_neg_equilibrium_test_passed << std::endl;
                //getchar();
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
                    //rates.val = rates_ref.val + drates;

                    // Fetch mole fractions
                    const auto& x_ref = properties_ref.moleFractions().val;
                    VectorConstRef xk_ref = x_ref(iks);
                    // TODO: loop for more then one kinetic species
                    bool kinetics_r_variation_check = true;
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
                //std::cout << "is_kin_rate_variation_test_passed = " << is_kin_rate_variation_test_passed << std::endl;
                //getchar();
                if(!is_kin_rate_variation_test_passed)
                    continue;

                result.timing.estimate_error_control += toc(ERROR_CONTROL_STEP);

                result.timing.estimate_search = toc(SEARCH_STEP);
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
                result.estimate.accepted = true;

                // Update the time
                t += dt;

                return;
            }
            else {
                inode_prev = inode;
                continue;
            }
        }

        result.estimate.accepted = false;

    }

    /// ------------------------------------------------------------------------------------------------------------ ///
    /// Solve functions
    /// ------------------------------------------------------------------------------------------------------------ ///

    auto solve(ChemicalState& state, double t, double dt) -> double
    {
        // Initialize result structure
        result = {};

        tic(SOLVE_STEP);

        //------------------------------------------------------------------------------------------
        // INITIALIZE OF THE VARIABLES DURING THE KINETICS SOLVE STEP
        //------------------------------------------------------------------------------------------
        tic(INITIALIZE_STEP);

        // Extract the composition of the kinetic species from the state
        const auto &n = state.speciesAmounts();
        ne = n(ies);
        nk = n(iks);

        // Assemble the vector benk = [be nk]
        benk.resize(Ee + Nk);
        benk.tail(Nk) = nk;

        result.timing.initialize = toc(INITIALIZE_STEP);

        //-----------------------------------------------------------------------------------------------
        // ESTIMATE ELEMENTS AND KINETICS' SPECIES STORED IN 'BENK' DURING THE SMART-KINETICS SOLVE STEP
        //-----------------------------------------------------------------------------------------------
        tic(ESTIMATE_STEP);

        // Perform a smart estimate for the chemical state
        if(options.method == SmartKineticStrategy::Clustering)
            estimate_clustering(state, t, dt);
        else if(options.method == SmartKineticStrategy::PriorityQueue)
            estimate_priority_based_acceptance_potential(state, t, dt);
        else if(options.method == SmartKineticStrategy::NearestNeighbour)
            estimate_nn_search_acceptance_based_lna(state, t, dt);

        result.timing.estimate = toc(ESTIMATE_STEP);

        //-----------------------------------------------------------------------------------------------
        // INTEGRATE ELEMENTS AND KINETICS' SPECIES STORED IN 'BENK' DURING THE SMART-KINETICS SOLVE STEP
        //-----------------------------------------------------------------------------------------------
        tic(LEARN_STEP);

        // Perform a learning step if the smart prediction is not satisfactory
        if(!result.estimate.accepted){
            if(options.method == SmartKineticStrategy::Clustering)
                learn_clustering(state, t, dt);
            else if(options.method == SmartKineticStrategy::PriorityQueue)
                learn_priority_based_acceptance_potential(state, t, dt);
            else if(options.method == SmartKineticStrategy::NearestNeighbour)
                learn_nn_search(state, t, dt);
        }

        // Extract the `be` and `nk` entries of the vector `benk`
        be = benk.head(Ee);
        nk = benk.tail(Nk);

        // Update the composition of the kinetic species
        state.setSpeciesAmounts(nk, iks);

        result.timing.learn = toc(LEARN_STEP);

        //------------------------------------------------------------------------------------------
        // EQUILIBRATE EQUILIBRIUM SPECIES STORED IN 'NE' DURING THE KINETICS SOLVE STEP
        //------------------------------------------------------------------------------------------

        tic(EQUILIBRATE_STEP);

        if(options.use_smart_equilibrium_solver){
            SmartEquilibriumResult res = {};
            res += smart_equilibrium.solve(state, T, P, be);
            result.smart_equilibrium += res;
        }
        else{
            EquilibriumResult res = {};
            res += equilibrium.solve(state, T, P, be);
            result.equilibrium += res;
        }

        result.timing.equilibrate = toc(EQUILIBRATE_STEP);

        result.timing.solve = toc(SOLVE_STEP);

        return t;

    }
    /// Solve the chemical kinetics problem from a given initial time to a final time.
    auto solve(ChemicalState& state, double t, double dt, VectorConstRef b) -> void
    {
        // Reset the result of the last smart equilibrium calculation
        result = {};

        tic(SOLVE_STEP);

        //------------------------------------------------------------------------------------------
        // INITIALIZE OF THE VARIABLES DURING THE KINETICS SOLVE STEP
        //------------------------------------------------------------------------------------------
        tic(INITIALIZE_STEP);

        // Extract the composition of the kinetic species from the state
        const auto &n = state.speciesAmounts();
        nk = n(iks);
        // Calculate the elements amounts related to the equilibrium species
        be = b - Ak * nk;

        // Assemble the vector benk = [be nk]
        benk.resize(Ee + Nk);
        benk.head(Ee) = be;
        benk.tail(Nk) = nk;

        // Initialise the chemical kinetics solver
        initialize(state, t, benk);

        result.timing.initialize=toc(INITIALIZE_STEP);

        //-----------------------------------------------------------------------------------------------
        // ESTIMATE ELEMENTS AND KINETICS' SPECIES STORED IN 'BENK' DURING THE SMART-KINETICS SOLVE STEP
        //-----------------------------------------------------------------------------------------------
        tic(ESTIMATE_STEP);

        // Perform a smart estimate for the chemical state
        if(options.method == SmartKineticStrategy::Clustering)
            estimate_clustering(state, t, dt);
        else if(options.method == SmartKineticStrategy::ClusteringExtended)
            estimate_clustering_extended(state, t, dt);
        else if(options.method == SmartKineticStrategy::PriorityQueue)
            estimate_priority_based_acceptance_potential(state, t, dt);
        else if(options.method == SmartKineticStrategy::PriorityQueuePrimary)
            estimate_priority_based_acceptance_primary_potential(state, t, dt);
        else if(options.method == SmartKineticStrategy::NearestNeighbour)
            estimate_nn_search_acceptance_based_lna(state, t, dt);

        result.timing.estimate = toc(ESTIMATE_STEP);

        //-----------------------------------------------------------------------------------------------
        // INTEGRATE ELEMENTS AND KINETICS' SPECIES STORED IN 'BENK' DURING THE SMART-KINETICS SOLVE STEP
        //-----------------------------------------------------------------------------------------------
        tic(LEARN_STEP);

        // Perform a learning step if the smart prediction is not satisfactory
        if(!result.estimate.accepted){
            if(options.method == SmartKineticStrategy::Clustering)
                learn_clustering(state, t, dt);
            else if(options.method == SmartKineticStrategy::ClusteringExtended)
                learn_clustering_extended(state, t, dt);
            else if(options.method == SmartKineticStrategy::PriorityQueue)
                learn_priority_based_acceptance_potential(state, t, dt);
            else if(options.method == SmartKineticStrategy::PriorityQueuePrimary)
                learn_priority_based_acceptance_primary_potential(state, t, dt);
            else if(options.method == SmartKineticStrategy::NearestNeighbour)
                learn_nn_search(state, t, dt);
        }

        // Extract the `be` and `nk` entries of the vector `benk`
        be = benk.head(Ee);
        nk = benk.tail(Nk);

        // Update the composition of the kinetic species
        state.setSpeciesAmounts(nk, iks);

        result.timing.learn = toc(LEARN_STEP);

        // ----------------------------------------------------------------------------------
        // Equilibrate equilibrium species
        // ----------------------------------------------------------------------------------
        tic(EQUILIBRATE_STEP);

        if(options.use_smart_equilibrium_solver){
            SmartEquilibriumResult res = {};
            res += smart_equilibrium.solve(state, T, P, be);
            result.smart_equilibrium += res;
        }
        else{
            EquilibriumResult res = {};
            res += equilibrium.solve(state, T, P, be);
            result.equilibrium += res;
        }

        result.timing.equilibrate = toc(EQUILIBRATE_STEP);

        result.timing.solve = toc(SOLVE_STEP);

    }

    auto function(ChemicalState& state, double t, VectorConstRef u, VectorRef res) -> int
    {
        // Extract the `be` and `nk` entries of the vector [be, nk]
        be = u.head(Ee);
        nk = u.tail(Nk);

        // Check for non-finite values in the vector `benk`
        for(Index i = 0; i < u.rows(); ++i)
            if(!std::isfinite(u[i]))
                return 1; // ensure the ode solver will reduce the time step

        // Update the composition of the kinetic species in the member `state`
        state.setSpeciesAmounts(nk, iks);

        // Solve the equilibrium problem using the elemental molar abundance `be`
        //if(options.use_smart_equilibrium_solver) // using smart equilibrium solver
        if(false) // NOTE: for now, using equilibrium solver
        {
            SmartEquilibriumResult result_eq = {};

            // smart_equilibrium_result
            timeit(result_eq += smart_equilibrium.solve(state, T, P, be), result.timing.learn_equilibration+=);

           result.smart_equilibrium += result_eq;

            // If smart calculation failed, use cold-start
            if(!result_eq.estimate.accepted && !result_eq.learning.gibbs_energy_minimization.optimum.succeeded)
            {
                //std::cout << "restart smart_equilibrium: " << res.learning.gibbs_energy_minimization.optimum.succeeded << std::endl;

                state.setSpeciesAmounts(0.0);
                timeit( result_eq = smart_equilibrium.solve(state, T, P, be), result.timing.learn_equilibration+=);

            }

            // Assert the smart equilibrium calculation did not fail
            Assert(result_eq.estimate.accepted || result_eq.learning.gibbs_energy_minimization.optimum.succeeded,
                   "Could not calculate the rates of the species.",
                   "The smart equilibrium calculation failed.");

        }
        else  // using conventional equilibrium solver
        {
            EquilibriumResult result_eq = {};

            timeit(result_eq += equilibrium.solve(state, T, P, be), result.timing.learn_equilibration +=);

            result.equilibrium += result_eq;

            // Check if the calculation failed, if so, use cold-start
            if(!result_eq.optimum.succeeded) {
                state.setSpeciesAmounts(0.0);
                timeit(result_eq = equilibrium.solve(state, T, P, be), result.timing.learn_equilibration +=);

            }
            if(!result_eq.optimum.succeeded) {
                return 1; // ensure the ode solver will reduce the time step
            }
            // Assert the equilibrium calculation did not fail
            Assert(result_eq.optimum.succeeded,
                   "Could not calculate the rates of the species.",
                   "The equilibrium calculation failed.");
        }

        // Update the chemical properties of the system
        timeit(properties = state.properties(), result.timing.learn_chemical_properties+=);

        // Calculate the kinetic rates of the reactions
        timeit(rates = reactions.rates(properties), result.timing.learn_reaction_rates+=);

        // Calculate the right-hand side function of the ODE
        res = A * rates.val;

        // Add the function contribution from the source rates
        if(source_fn)
        {
            // Evaluate the source function
            q = source_fn(properties);

            // Add the contribution of the source rates
            res += B * q.val;
        }

        return 0;
    }

    auto jacobian(ChemicalState& state, double t, VectorConstRef u, MatrixRef res) -> int
    {
        // Calculate the sensitivity of the equilibrium state
        //timeit( sensitivity = options.use_smart_equilibrium_solver ? smart_equilibrium.sensitivity() : equilibrium.sensitivity(),
        //        result.timing.learn_sensitivity+=);
        sensitivity = equilibrium.sensitivity();

        // Extract the columns of the kinetic rates derivatives w.r.t. the equilibrium and kinetic species
        drdne = cols(rates.ddn, ies);
        drdnk = cols(rates.ddn, iks);

        // Calculate the derivatives of `r` w.r.t. `be` using the equilibrium sensitivity
        drdbe = drdne * sensitivity.dndb(ies, iee);

        // Assemble the partial derivatives of the reaction rates `r` w.r.t. to `u = [be nk]`
        drdu << drdbe, drdnk;

        // Calculate the Jacobian matrix of the ODE function
        res = A * drdu;

        // Add the Jacobian contribution from the source rates
        if(source_fn)
        {
            // Extract the columns of the source rates derivatives w.r.t. the equilibrium and kinetic species
            dqdne = cols(q.ddn, ies);
            dqdnk = cols(q.ddn, iks);

            // Calculate the derivatives of `q` w.r.t. `be` using the equilibrium sensitivity
            dqdbe = dqdne * sensitivity.dndb(ies, iee);

            // Assemble the partial derivatives of the source rates `q` w.r.t. to `u = [be nk]`
            dqdu << dqdbe, dqdnk;

            // Add the contribution of the source rates
            res += B * dqdu;
        }

        return 0;
    }

    auto outputSmartMethodInfo() const -> void
    {
        std::cout << "***********************************************************************************" << std::endl;
        std::cout << "Clusters ordered by order of their creation" << std::endl;
        std::cout << "***********************************************************************************" << std::endl;

        Index i = 0;
        if(options.method == SmartKineticStrategy::Clustering)
            for(auto cluster : database.clusters)
            {
                std::cout << "CLUSTER #" << i << std::endl;
                std::cout << "  RANK OF CLUSTER: " << database.priority.priorities()[i] << std::endl;
                std::cout << "  PRIMARY SPECIES: ";
                for(auto j : database.clusters[i].iprimary)
                    std::cout << system.species(j).name() << " ";
                std::cout << std::endl;
                std::cout << "  NUMBER OF RECORDS: " << database.clusters[i].records.size() << std::endl;
                //std::cout << "  RANK OF RECORDS: ";
                //for(auto j : database.clusters[i].priority.order())
                //    std::cout << database.clusters[i].priority.priorities()[j] << " ";
                std::cout << std::endl;
                std::cout << std::endl;
                i++;
            }
        else if (options.method == SmartKineticStrategy::ClusteringExtended)
            for(auto cluster : database_extended.clusters)
            {
                std::cout << "CLUSTER #" << i << std::endl;
                std::cout << "  RANK OF CLUSTER: " << database_extended.priority.priorities()[i] << std::endl;
                std::cout << "  PRIMARY SPECIES: ";
                for(auto j : database_extended.clusters[i].iprimary_ikin)
                    std::cout << system.species(j).name() << " ";
                std::cout << std::endl;
                std::cout << "  NUMBER OF RECORDS: " << database_extended.clusters[i].records.size() << std::endl;
                //std::cout << "  RANK OF RECORDS: ";
                //for(auto j : database.clusters[i].priority.order())
                //    std::cout << database.clusters[i].priority.priorities()[j] << " ";
                std::cout << std::endl;
                std::cout << std::endl;
                i++;
            }
        std::cout << "***********************************************************************************" << std::endl;

        // Output content of the smart equilibrium solver if it is used
        if(options.use_smart_equilibrium_solver)
            smart_equilibrium.outputInfo();

        /*
        for(auto i : database.priority.order())
        {
            std::cout << "CLUSTER #" << i << std::endl;
            std::cout << "  RANK OF CLUSTER: " << database.priority.priorities()[i] << std::endl;
            std::cout << "  PRIMARY SPECIES: ";
            for(auto j : database.clusters[i].iprimary)
                std::cout << system.species(j).name() << " ";
            std::cout << std::endl;
            std::cout << "  NUMBER OF RECORDS: " << database.clusters[i].records.size() << std::endl;
            std::cout << "  RANK OF RECORDS: ";
            for(auto j : database.clusters[i].priority.order())
                std::cout << database.clusters[i].priority.priorities()[j] << " ";
            std::cout << std::endl;
            std::cout << "  NEXT CLUSTER: " << database.connectivity.order(i)[1] << std::endl;
            std::cout << "  NEXT CLUSTER PRIMARY SPECIES: ";
            for(auto j : database.clusters[database.connectivity.order(i)[1]].iprimary)
                std::cout << system.species(j).name() << " ";
            std::cout << std::endl;
            std::cout << std::endl;
        }
        */
    }
};

SmartKineticSolver::SmartKineticSolver()
{
    RuntimeError("Cannot proceed with SmartKineticSolver().",
                 "SmartKineticSolver() constructor is deprecated. "
                 "Use constructor SmartKineticSolver(const ReactionSystem&, const Partition&) instead.");
}

SmartKineticSolver::SmartKineticSolver(const ReactionSystem& reactions)
{
    RuntimeError("Cannot proceed with SmartKineticSolver(const ReactionSystem&).",
                 "SmartKineticSolver(const ReactionSystem&) constructor is deprecated. "
                 "Use constructor SmartKineticSolver(const ReactionSystem&, const Partition&) instead.");
}


SmartKineticSolver::SmartKineticSolver(const ReactionSystem& reactions, const Partition& partition)
: pimpl(new Impl(reactions, partition))
{}


SmartKineticSolver::~SmartKineticSolver()
{
}

auto SmartKineticSolver::operator=(SmartKineticSolver other) -> SmartKineticSolver&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto SmartKineticSolver::setOptions(const SmartKineticOptions& options) -> void
{
    pimpl->setOptions(options);
}

auto SmartKineticSolver::setPartition(const Partition& partition) -> void
{
    RuntimeError("Cannot proceed with SmartKineticSolver::setPartition.",
                 "SmartKineticSolver::setPartition is deprecated. "
                 "Use constructor SmartKinetic"
                 "Solver(const Partition&) instead.");
}

auto SmartKineticSolver::addSource(const ChemicalState& state, double volumerate, const std::string& units) -> void
{
    pimpl->addSource(state, volumerate, units);
}

auto SmartKineticSolver::addPhaseSink(const std::string& phase, double volumerate, const std::string& units) -> void
{
    pimpl->addPhaseSink(phase, volumerate, units);
}

auto SmartKineticSolver::addFluidSink(double volumerate, const std::string& units) -> void
{
    pimpl->addFluidSink(volumerate, units);
}

auto SmartKineticSolver::addSolidSink(double volumerate, const std::string& units) -> void
{
    pimpl->addSolidSink(volumerate, units);
}

auto SmartKineticSolver::initialize(ChemicalState& state, double tstart) -> void
{
    pimpl->initialize(state, tstart);
}

auto SmartKineticSolver::initialize(ChemicalState& state, double tstart, VectorConstRef benk) -> void
{
    pimpl->initialize(state, tstart, benk);
}

auto SmartKineticSolver::result() const -> const SmartKineticResult&
{
    return pimpl->result;
}

auto SmartKineticSolver::properties() const -> const ChemicalProperties&
{
    return pimpl->properties;
}

auto SmartKineticSolver::solve(ChemicalState& state, double t, double dt) -> double
{
    return pimpl->solve(state, t, dt);
}

auto SmartKineticSolver::solve(ChemicalState& state, double t, double dt, VectorConstRef b) -> void
{
    pimpl->solve(state, t, dt, b);
}

auto SmartKineticSolver::outputSmartMethodInfo() const -> void
{
    return pimpl->outputSmartMethodInfo();
}

} // namespace Reaktoro
