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

#include "SmartKineticSolverBase.hpp"

namespace Reaktoro {

SmartKineticSolverBase::SmartKineticSolverBase(const ReactionSystem& reactions, const Partition& partition)
: _reactions(reactions), system(partition.system()), _partition(partition),
  equilibrium(partition), smart_equilibrium(partition)
{
    // Set the indices of the equilibrium, kinetic, and inert species
    ies = _partition.indicesEquilibriumSpecies();
    iks = _partition.indicesKineticSpecies();

    // Set the indices of the equilibrium and kinetic elements
    iee = _partition.indicesEquilibriumElements();
    ike = _partition.indicesKineticElements();

    // Set the number of equilibrium, kinetic, and inert species
    Ne = ies.size();
    Nk = iks.size();
    N = system.numSpecies();

    // Set the number of equilibrium and kinetic elements
    Ee = iee.size();
    E = system.numElements();

    // Initialise the formula matrix of the equilibrium partition
    Ae = _partition.formulaMatrixEquilibriumPartition();

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
    Se = cols(_reactions.stoichiometricMatrix(), ies);
    Sk = cols(_reactions.stoichiometricMatrix(), iks);

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
    drdu.resize(_reactions.numReactions(), Ee + Nk);

    // Allocate memory for the partial derivatives of the source rates `q` w.r.t. to `u = [be nk]`
    dqdu.resize(N, Ee + Nk);

    // Allocate the memory for the sensitivity matrix
    benk_S.resize(Ee + Nk, Ee + Nk);
}

SmartKineticSolverBase::~SmartKineticSolverBase()
{
}

auto SmartKineticSolverBase::setOptions(const SmartKineticOptions& _options) -> void
{
    // Initialise the options of the kinetic solver
    this->options = _options;

    // Set options for the equilibrium calculations
    equilibrium.setOptions(options.equilibrium);
    smart_equilibrium.setOptions(options.smart_equilibrium);

}

auto SmartKineticSolverBase::setPartition(const Partition& _partition) -> void
{
    RuntimeError("Cannot proceed with SmartKineticSolverBase::setPartition.",
                 "SmartKineticSolverBase::setPartition is deprecated. "
                 "Use constructor SmartKinetic"
                 "Solver(const Partition&) instead.")
}

auto SmartKineticSolverBase::addSource(ChemicalState state, double volumerate, const std::string& units) -> void
{
    const Index num_species = system.numSpecies();
    const double volume = units::convert(volumerate, units, "m3/s");
    state.scaleVolume(volume);
    const Vector _n = state.speciesAmounts();
    auto old_source_fn = source_fn;

    source_fn = [=](const ChemicalProperties& properties)
    {
        ChemicalVector _q(num_species);
        _q.val = _n;
        if(old_source_fn)
            _q += old_source_fn(properties);
        return _q;
    };
}

auto SmartKineticSolverBase::addPhaseSink(const std::string& phase, double volumerate, const std::string& units) -> void
{
    const double volume = units::convert(volumerate, units, "m3/s");
    const Index iphase = system.indexPhaseWithError(phase);
    const Index ifirst = system.indexFirstSpeciesInPhase(iphase);
    const Index size = system.numSpeciesInPhase(iphase);
    auto old_source_fn = source_fn;
    ChemicalScalar phasevolume;
    ChemicalVector _q(size);

    source_fn = [=](const ChemicalProperties& properties) mutable
    {
        const auto _n = properties.composition();
        const auto np = rows(_n, ifirst, size);
        auto qp = rows(_q, ifirst, size);
        phasevolume = properties.phaseVolumes()[iphase];
        qp = -volume*np/phasevolume;
        if(old_source_fn)
            _q += old_source_fn(properties);
        return _q;
    };
}

auto SmartKineticSolverBase::addFluidSink(double volumerate, const std::string& units) -> void
{
    const double volume = units::convert(volumerate, units, "m3/s");
    const Indices& isolid_species = _partition.indicesSolidSpecies();
    auto old_source_fn = source_fn;
    ChemicalScalar fluidvolume;
    ChemicalVector _q;

    source_fn = [=](const ChemicalProperties& properties) mutable
    {
        const auto _n = properties.composition();
        fluidvolume = properties.fluidVolume();
        _q = -volume*_n/fluidvolume;
        rows(_q, isolid_species).fill(0.0);
        if(old_source_fn)
            _q += old_source_fn(properties);
        return _q;
    };
}

auto SmartKineticSolverBase::addSolidSink(double volumerate, const std::string& units) -> void
{
    const double volume = units::convert(volumerate, units, "m3/s");
    const Indices& ifluid_species = _partition.indicesFluidSpecies();
    auto old_source_fn = source_fn;
    ChemicalScalar solidvolume;
    ChemicalVector _q;

    source_fn = [=](const ChemicalProperties& properties) mutable
    {
        const auto _n = properties.composition();
        solidvolume = properties.solidVolume();
        _q = -volume*_n/solidvolume;
        rows(_q, ifluid_species).fill(0.0);
        if(old_source_fn)
            _q += old_source_fn(properties);
        return _q;
    };
}

auto SmartKineticSolverBase::initialize( ChemicalState& state, double tstart) -> void {

    // Extract the composition of the equilibrium and kinetic species
    const auto &_n = state.speciesAmounts();
    nk = _n(iks);
    ne = _n(ies);

    // Assemble the vector benk = [be nk]
    benk.resize(Ee + Nk);
    benk.head(Ee) = Ae * ne;
    benk.tail(Nk) = nk;

    // Initialize
    initialize(state, tstart, benk);

}

auto SmartKineticSolverBase::initialize(ChemicalState& state, double tstart, VectorConstRef _benk) -> void
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
    ode.initialize(tstart, _benk);

    // Set the options of the equilibrium solver
    smart_equilibrium.setOptions(options.learning.smart_equilibrium);
    equilibrium.setOptions(options.learning.equilibrium);

}

auto SmartKineticSolverBase::function(ChemicalState& state, double t, VectorConstRef u, VectorRef res) -> int
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
        timeit(result_eq += smart_equilibrium.solve(state, T, P, be), _result.timing.learn_equilibration+=)

        _result.smart_equilibrium += result_eq;

        // If smart calculation failed, use cold-start
        if(!result_eq.estimate.accepted && !result_eq.learning.gibbs_energy_minimization.optimum.succeeded)
        {
            state.setSpeciesAmounts(0.0);
            timeit( result_eq = smart_equilibrium.solve(state, T, P, be), _result.timing.learn_equilibration+=)

        }

        // Assert the smart equilibrium calculation did not fail
        Assert(result_eq.estimate.accepted || result_eq.learning.gibbs_energy_minimization.optimum.succeeded,
               "Could not calculate the rates of the species.",
               "The smart equilibrium calculation failed.")

    }
    else  // using conventional equilibrium solver
    {
        EquilibriumResult result_eq = {};

        timeit(result_eq += equilibrium.solve(state, T, P, be), _result.timing.learn_equilibration +=)

        _result.equilibrium += result_eq;

        // Check if the calculation failed, if so, use cold-start
        if(!result_eq.optimum.succeeded) {
            state.setSpeciesAmounts(0.0);
            timeit(result_eq = equilibrium.solve(state, T, P, be), _result.timing.learn_equilibration +=)

        }
        if(!result_eq.optimum.succeeded) {
            return 1; // ensure the ode solver will reduce the time step
        }
        // Assert the equilibrium calculation did not fail
        Assert(result_eq.optimum.succeeded,
               "Could not calculate the rates of the species.",
               "The equilibrium calculation failed.")
    }

    // Update the chemical properties of the system
    timeit(_properties = state.properties(), _result.timing.learn_chemical_properties+=)

    // Calculate the kinetic rates of the reactions
    timeit(rates = _reactions.rates(_properties), _result.timing.learn_reaction_rates+=)

    // Calculate the right-hand side function of the ODE
    res = A * rates.val;

    // Add the function contribution from the source rates
    if(source_fn)
    {
        // Evaluate the source function
        q = source_fn(_properties);

        // Add the contribution of the source rates
        res += B * q.val;
    }

    return 0;
}

auto SmartKineticSolverBase::jacobian(ChemicalState& state, double t, VectorConstRef u, MatrixRef res) -> int
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

auto SmartKineticSolverBase::result() const -> const SmartKineticResult&
{
    return _result;
}

auto SmartKineticSolverBase::properties() const -> const ChemicalProperties&
{
    return _properties;
}

auto SmartKineticSolverBase::reactions() const -> const ReactionSystem&
{
    return _reactions;
}

auto SmartKineticSolverBase::partition() const -> const Partition&
{
    return _partition;
}

auto SmartKineticSolverBase::solve(ChemicalState& state, double t, double dt) -> double
{
    // Initialize result structure
    _result = {};

    tic(SOLVE_STEP)

    //------------------------------------------------------------------------------------------
    // INITIALIZE OF THE VARIABLES DURING THE KINETICS SOLVE STEP
    //------------------------------------------------------------------------------------------
    tic(INITIALIZE_STEP)

    // Extract the composition of the kinetic species from the state
    const auto &_n = state.speciesAmounts();
    ne = _n(ies);
    nk = _n(iks);

    // Assemble the vector benk = [be nk]
    //benk.head(Ee) = be;
    benk.head(Ee) = Ae * ne;
    benk.tail(Nk) = nk;

    _result.timing.initialize = toc(INITIALIZE_STEP);

    //-----------------------------------------------------------------------------------------------
    // ESTIMATE ELEMENTS AND KINETICS' SPECIES STORED IN 'BENK' DURING THE SMART-KINETICS SOLVE STEP
    //-----------------------------------------------------------------------------------------------
    tic(ESTIMATE_STEP)

    // Perform a smart estimate for the chemical state
    estimate(state, t, dt);

    _result.timing.estimate = toc(ESTIMATE_STEP);

    //-----------------------------------------------------------------------------------------------
    // INTEGRATE ELEMENTS AND KINETICS' SPECIES STORED IN 'BENK' DURING THE SMART-KINETICS SOLVE STEP
    //-----------------------------------------------------------------------------------------------
    tic(LEARN_STEP)

    // Perform a learning step if the smart prediction is not satisfactory
    if(!_result.estimate.accepted){
        learn(state, t, dt);
    }

    // Extract the `be` and `nk` entries of the vector `benk`
    be = benk.head(Ee);
    nk = benk.tail(Nk);

    // Update the composition of the kinetic species
    state.setSpeciesAmounts(nk, iks);

    _result.timing.learn = toc(LEARN_STEP);

    //------------------------------------------------------------------------------------------
    // EQUILIBRATE EQUILIBRIUM SPECIES STORED IN 'NE' DURING THE KINETICS SOLVE STEP
    //------------------------------------------------------------------------------------------

    tic(EQUILIBRATE_STEP)

    if(options.use_smart_equilibrium_solver){
        SmartEquilibriumResult res = {};
        res += smart_equilibrium.solve(state, T, P, be);
        _result.smart_equilibrium += res;
    }
    else{
        EquilibriumResult res = {};
        res += equilibrium.solve(state, T, P, be);
        _result.equilibrium += res;
    }

    _result.timing.equilibrate = toc(EQUILIBRATE_STEP);

    _result.timing.solve = toc(SOLVE_STEP);

    return t;

}

auto SmartKineticSolverBase::solve(ChemicalState& state, double t, double dt, VectorConstRef b) -> void
{
    // Reset the result of the last smart equilibrium calculation
    _result = {};

    tic(SOLVE_STEP)

    //------------------------------------------------------------------------------------------
    // INITIALIZE OF THE VARIABLES DURING THE KINETICS SOLVE STEP
    //------------------------------------------------------------------------------------------
    tic(INITIALIZE_STEP)

    // Extract the composition of the kinetic species from the state
    const auto &_n = state.speciesAmounts();
    nk = _n(iks);
    // Calculate the elements amounts related to the equilibrium species
    be = b - Ak * nk;

    // Assemble the vector benk = [be nk]
    benk.resize(Ee + Nk);
    benk.head(Ee) = be;
    benk.tail(Nk) = nk;

    // Initialise the chemical kinetics solver
    initialize(state, t, benk);

    _result.timing.initialize = toc(INITIALIZE_STEP);

    //-----------------------------------------------------------------------------------------------
    // ESTIMATE ELEMENTS AND KINETICS' SPECIES STORED IN 'BENK' DURING THE SMART-KINETICS SOLVE STEP
    //-----------------------------------------------------------------------------------------------
    tic(ESTIMATE_STEP)

    // Perform a smart estimate for the chemical state
    estimate(state, t, dt);

    _result.timing.estimate = toc(ESTIMATE_STEP);

    //-----------------------------------------------------------------------------------------------
    // INTEGRATE ELEMENTS AND KINETICS' SPECIES STORED IN 'BENK' DURING THE SMART-KINETICS SOLVE STEP
    //-----------------------------------------------------------------------------------------------
    tic(LEARN_STEP)

    // Perform a learning step if the smart prediction is not satisfactory
    learn(state, t, dt);

    // Extract the `be` and `nk` entries of the vector `benk`
    be = benk.head(Ee);
    nk = benk.tail(Nk);

    // Update the composition of the kinetic species
    state.setSpeciesAmounts(nk, iks);

    _result.timing.learn = toc(LEARN_STEP);

    // ----------------------------------------------------------------------------------
    // Equilibrate equilibrium species
    // ----------------------------------------------------------------------------------
    tic(EQUILIBRATE_STEP)

    if(options.use_smart_equilibrium_solver){
        SmartEquilibriumResult res = {};
        res += smart_equilibrium.solve(state, T, P, be);
        _result.smart_equilibrium += res;
    }
    else{
        EquilibriumResult res = {};
        res += equilibrium.solve(state, T, P, be);
        _result.equilibrium += res;
    }

    _result.timing.equilibrate = toc(EQUILIBRATE_STEP);

    _result.timing.solve = toc(SOLVE_STEP);

}


} // namespace Reaktoro
