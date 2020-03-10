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

#include "EquilibriumSolver.hpp"

// Reaktoro includes
#include <Reaktoro/Common/ChemicalVector.hpp>
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/ConvertUtils.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Profiling.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Connectivity.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Core/ThermoProperties.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/Math/MathUtils.hpp>
#include <Reaktoro/Optimization/OptimumOptions.hpp>
#include <Reaktoro/Optimization/OptimumProblem.hpp>
#include <Reaktoro/Optimization/OptimumResult.hpp>
#include <Reaktoro/Optimization/OptimumSolver.hpp>
#include <Reaktoro/Optimization/OptimumSolverRefiner.hpp>
#include <Reaktoro/Optimization/OptimumState.hpp>

namespace Reaktoro {

struct EquilibriumSolver::Impl
{
    /// The chemical system instance
    ChemicalSystem system;

    /// The partition of the chemical system
    Partition partition;

    /// The options of the equilibrium solver
    EquilibriumOptions options;

    /// The result of the last equilibrium calculation.
    EquilibriumResult result;

    /// The solver for the optimisation calculations
    OptimumSolver solver;

    /// The chemical properties of the chemical system
    ChemicalProperties properties;

    /// The sensitivity derivatives of the equilibrium state
    EquilibriumSensitivity sensitivities, sensitivities_e;

    /// The molar amounts of the species
    Vector n;

    /// The dual potentials of the elements
    Vector y;

    /// The dual potentials of the species
    Vector z;

    /// The molar amounts of the elements in the equilibrium partition
    Vector be;

    /// The standard chemical potentials of the species normalized by RT
    ThermoVector u0;

    /// The chemical potentials of the species normalized by RT
    ChemicalVector u;

    /// The chemical potentials of the equilibrium species normalized by RT
    ChemicalVector ue;

    /// The chemical potentials of the inert species normalized by RT
    Vector ui;

    /// The mole fractions of the equilibrium species
    ChemicalVector xe;

    /// The optimisation problem
    OptimumProblem optimum_problem;

    /// The state of the optimisation calculation
    OptimumState optimum_state;

    /// The options for the optimisation calculation
    OptimumOptions optimum_options;

    /// The indices of the species in the equilibrium partition
    Indices ies;

    /// The indices of the elements in the equilibrium partition
    Indices iee;

    /// The indices of the inert species (i.e., the species in disequilibrium)
    Indices iis;

    /// The indices of the fluid species (i.e., the species in disequilibrium)
    Indices ifs;

    /// The indices of the solid species (i.e., the species in disequilibrium)
    Indices iss;

    /// The number of species and elements in the system
    unsigned N, E;

    /// The number of species and elements in the equilibrium partition
    unsigned Ne, Ee;

    /// The formula matrix of the species in the system
    Matrix A;

    /// The formula matrix of the species in the equilibrium partition
    Matrix Ae;

    /// The formula matrix of the inert species
    Matrix Ai;

    int skipped = 0;

    /// Construct a default Impl instance
    Impl()
    {}

    /// Construct a Impl instance with given Partition
    Impl(const Partition& partition)
    : system(partition.system()), properties(partition.system())
    {
        // Initialize the formula matrix
        A = system.formulaMatrix();

        // Initialize the number of species and elements in the system
        N = system.numSpecies();
        E = system.numElements();

        // Initialize the number of species and elements in the equilibrium partition
        Ne = partition.numEquilibriumSpecies();
        Ee = partition.numEquilibriumElements();

        // Initialize the formula matrix of the equilibrium species
        Ae = partition.formulaMatrixEquilibriumPartition();

        // Initialize the indices of the equilibrium species and elements
        ies = partition.indicesEquilibriumSpecies();
        iee = partition.indicesEquilibriumElements();
        ifs = partition.indicesEquilibriumFluidSpecies();
        iss = partition.indicesEquilibriumSolidSpecies();

        // Initialize the indices of the inert species
        iis.clear();
        iis.reserve(partition.numInertSpecies() + partition.numKineticSpecies());
        iis.insert(iis.end(), partition.indicesInertSpecies().begin(), partition.indicesInertSpecies().end());
        iis.insert(iis.end(), partition.indicesKineticSpecies().begin(), partition.indicesKineticSpecies().end());

        // Initialize the formula matrix of the inert species
        Ai = cols(A, iis);

        // Initialize the sensitivity matrices with respect to all species in the system
        sensitivities.dndT = zeros(N);
        sensitivities.dndP = zeros(N);
        sensitivities.dndb = zeros(N, E);
        sensitivities.dydT = zeros(E);
        sensitivities.dydP = zeros(E);
        sensitivities.dydb = zeros(E, E);
        sensitivities.dzdT = zeros(N);
        sensitivities.dzdP = zeros(N);
        sensitivities.dzdb = zeros(N, E);

        // Initialize the sensitivity matrices with respect to equilibrium species only
        sensitivities_e.dndT = zeros(Ne);
        sensitivities_e.dndP = zeros(Ne);
        sensitivities_e.dndb = zeros(Ne, Ee);
        sensitivities_e.dydT = zeros(Ee);
        sensitivities_e.dydP = zeros(Ee);
        sensitivities_e.dydb = zeros(Ee, Ee);
        sensitivities_e.dzdT = zeros(Ne);
        sensitivities_e.dzdP = zeros(Ne);
        sensitivities_e.dzdb = zeros(Ne, Ee);
    }

    /// Update the OptimumOptions instance with given EquilibriumOptions instance
    auto updateOptimumOptions() -> void
    {
        // Initialize the options for the optimisation calculation
        optimum_options = options.optimum;

        // Set the parameters of the optimisation algorithms that control how small can be the amount of a species
        optimum_options.ipaction.mu = options.epsilon;
        optimum_options.ipnewton.mu = options.epsilon;
        optimum_options.ipopt.mu.push_back(options.epsilon);
        optimum_options.ipactive.epsilon = options.epsilon;

        // Initialize the names of the primal and dual variables
        if(options.optimum.output.active)
        {
            // Use `n` instead of `x` to name the variables
            optimum_options.output.xprefix = "n";

            // Define some auxiliary references to the variables names
            auto& xnames = optimum_options.output.xnames;
            auto& ynames = optimum_options.output.ynames;
            auto& znames = optimum_options.output.znames;

            // Initialize the names of the primal variables `n`
            for(Index i : ies)
                xnames.push_back(system.species(i).name());

            // Initialize the names of the dual variables `y`
            for(Index i : iee)
                ynames.push_back(system.element(i).name());

            // Initialize the names of the dual variables `z`
            znames = xnames;
        }
    }

    /// Update the OptimumProblem instance with given EquilibriumProblem and ChemicalState instances
    auto updateOptimumProblem(const ChemicalState& state) -> void
    {
        // The temperature and pressure of the equilibrium calculation
        const auto T  = state.temperature();
        const auto P  = state.pressure();
        const auto RT = universalGasConstant*T;

        // Set the molar amounts of the species
        n = state.speciesAmounts();

        // Update the standard thermodynamic properties of the chemical system
        timeit( properties.update(T, P),
            result.timing.standard_thermodynamic_properties+= );

        // Update the normalized standard Gibbs energies of the species
        u0 = properties.standardPartialMolarGibbsEnergies()/RT;

        // The result of the objective evaluation
        ObjectiveResult obj;

        // The Gibbs energy function to be minimized
        optimum_problem.objective = [=](VectorConstRef ne) mutable
        {
            // Set the molar amounts of the species
            n(ies) = ne;

            // Update the chemical properties of the chemical system
            timeit( properties.update(n),
                result.timing.chemical_properties+= );

            // Set the scaled chemical potentials of the species
            u = u0 + properties.lnActivities();

            // Set the scaled chemical potentials of the equilibrium species
            ue = rows(u, ies, ies);

            // Set the mole fractions of the equilibrium species
            xe = rows(properties.moleFractions(), ies, ies);

            // Set the objective result
            obj.val = dot(ne, ue.val);
            obj.grad = ue.val;

            // Set the Hessian of the objective function
            switch(options.hessian)
            {
            case GibbsHessian::Exact:
                obj.hessian.mode = Hessian::Dense;
                obj.hessian.dense = ue.ddn;
                break;
            case GibbsHessian::ExactDiagonal:
                obj.hessian.mode = Hessian::Diagonal;
                obj.hessian.diagonal = diagonal(ue.ddn);
                break;
            case GibbsHessian::Approximation:
                obj.hessian.mode = Hessian::Dense;
                obj.hessian.dense = diag(inv(xe.val)) * xe.ddn;
                break;
            case GibbsHessian::ApproximationDiagonal:
                obj.hessian.mode = Hessian::Diagonal;
                obj.hessian.diagonal = diagonal(xe.ddn)/xe.val;
                break;
            }

            return obj;
        };

        optimum_problem.c.resize(0);
        optimum_problem.n = Ne;
        optimum_problem.A = Ae;
        optimum_problem.b = be;
        optimum_problem.l.setConstant(Ne, options.epsilon);
    }

    /// Initialize the optimum state from a chemical state
    auto updateOptimumState(const ChemicalState& state) -> void
    {
        // The temperature and the RT factor
        const double T  = state.temperature();
        const double RT = universalGasConstant*T;

        // Set the molar amounts of the species
        n = state.speciesAmounts();

        // Set the normalized dual potentials of the elements
        y = state.equilibrium().elementChemicalPotentials()/RT;

        // Set the normalized dual potentials of the species
        z = state.equilibrium().speciesStabilities()/RT;

        // Initialize the optimum state
        optimum_state.x = n(ies);
        optimum_state.y = y(iee);
        optimum_state.z = z(ies);
    }

    /// Initialize the chemical state from a optimum state
    auto updateChemicalState(ChemicalState& state) -> void
    {
        // The temperature and the RT factor
        const double T  = state.temperature();
        const double RT = universalGasConstant*T;

        // Update the molar amounts of the equilibrium species
        n(ies) = optimum_state.x;

        // Update the normalized chemical potentials of the inert species
        ui = u.val(iis);

        // Update the normalized dual potentials of the elements
        y = zeros(E); y(iee) = optimum_state.y;

        // Update the normalized dual potentials of the equilibrium and inert species
        z(ies) = optimum_state.z;
        z(iis) = ui - tr(Ai) * y;

        // Scale the normalized dual potentials of elements and species to units of J/mol
        y *= RT;
        z *= RT;

        // Update the chemical state
        state.setSpeciesAmounts(n);
        state.equilibrium().setElementChemicalPotentials(y);
        state.equilibrium().setSpeciesStabilities(z);
        state.equilibrium().setIndicesStrictlyUnstableElements(optimum_state.itrivial_constraints);
        state.equilibrium().setIndicesStrictlyUnstableSpecies(optimum_state.itrivial_variables);
    }

    /// Find a feasible approximation for an equilibrium problem with all elements present on chemical system.
    auto approximate_with_all_element_amounts(ChemicalState& state, double T, double P, VectorConstRef b) -> EquilibriumResult
    {
        Vector be = b(iee);
        return approximate(state, T, P, be);
    }

    /// Find a feasible approximation for an equilibrium problem.
    auto approximate(ChemicalState& state, double T, double P, Vector be) -> EquilibriumResult
    {
        // Check the dimension of the vector `be`
        Assert(unsigned(be.rows()) == Ee,
            "Cannot proceed with method EquilibriumSolver::approximate.",
            "The dimension of the given vector of molar amounts of the "
            "elements does not match the number of elements in the "
            "equilibrium partition.");

        // Set temperature and pressure of the chemical state
        state.setTemperature(T);
        state.setPressure(P);

        // Auxiliary variables
        const double RT = universalGasConstant*T;
        const double inf = std::numeric_limits<double>::infinity();

        // Update the internal state of n, y, z
        n = state.speciesAmounts();
        y = state.equilibrium().elementChemicalPotentials();
        z = state.equilibrium().speciesStabilities();

        // Update the standard thermodynamic properties of the system
        timeit( properties.update(T, P),
                result.timing.standard_thermodynamic_properties+= );

        // Get the standard Gibbs energies of the equilibrium species
        const Vector ge0 = properties.standardPartialMolarGibbsEnergies().val(ies);

        // Get the ln activity constants of the equilibrium species
        const Vector ln_ce = properties.lnActivityConstants().val(ies);

        // Define the optimisation problem
        OptimumProblem optimum_problem;
        optimum_problem.n = Ne;
        optimum_problem.c = ge0/RT + ln_ce;
        optimum_problem.A = Ae;
        optimum_problem.b = be;
        optimum_problem.l = zeros(Ne);
        optimum_problem.u = ones(Ne) * inf;

        // Initialize the optimum state
        OptimumState optimum_state;

        // The result of the linear programming calculation
        EquilibriumResult result;

        // Update the optimum options
        updateOptimumOptions();

        // Set the method for the optimisation calculation
        solver.setMethod(options.method);

        // Solve the linear programming problem
        result.optimum = solver.solve(optimum_problem, optimum_state, optimum_options);

        // Update the molar amounts of the equilibrium species
        n(ies) = optimum_state.x;

        // Update the dual potentials of the species and elements (in units of J/mol)
        z = zeros(N); z(ies) = optimum_state.z * RT;
        y = zeros(E); y(iee) = optimum_state.y * RT;

        // Update the chemical state
        state.setSpeciesAmounts(n);
        state.equilibrium().setElementChemicalPotentials(y);
        state.equilibrium().setSpeciesStabilities(z);

        return result;
    }

    /// Find an initial guess for an equilibrium problem.
    auto initialguess(ChemicalState& state, double T, double P, Vector be) -> EquilibriumResult
    {
        // Solve the linear programming problem to obtain an approximation
        auto result = approximate(state, T, P, be);

        // Check the approximate calculation was successful
        Assert(result.optimum.succeeded,
            "Cannot proceed with the equilibrium calculation.",
            "The calculation of initial guess failed, most "
            "probably because no feasible solution exists.");

        // Preserve the values of n that are greater than z. For all others, set n to sqrt(epsilon)
        n = (n.array() > z.array()).select(n, std::sqrt(options.epsilon));

        // Set z to zero, which will later be set to epsilon / n.
        z.fill(0.0);

        // Update the chemical state
        state.setSpeciesAmounts(n);
        state.equilibrium().setSpeciesStabilities(z);

        return result;
    }

    /// Return true if cold-start is needed.
    auto coldstart(const ChemicalState& state) -> bool
    {
        // Check if all equilibrium species have zero amounts
        bool zero = true;
        for(Index i : ies)
            if(state.speciesAmount(i) > 0)
                { zero = false; break; }
        return zero || !options.warmstart;
    }

    /// Solve the equilibrium problem with given initial state
    auto solve(ChemicalState& state) -> EquilibriumResult
    {
        const auto T = state.temperature();
        const auto P = state.pressure();
        be = state.elementAmountsInSpecies(ies)(iee);
        return solve(state, T, P, be);
    }

    /// Solve the equilibrium problem with given problem definition
    auto solve(ChemicalState& state, const EquilibriumProblem& problem) -> EquilibriumResult
    {
        Assert(problem.partition() == partition,
            "Cannot proceed with method EquilibriumSolver::solve(ChemicalState&, const EquilibriumProblem&).",
            "This EquilibriumSolver and the given EquilibriumProblem objects were not initialized with the same Partition object.");
        const auto T = problem.temperature();
        const auto P = problem.pressure();
        be = problem.elementAmounts()(iee);
        return solve(state, T, P, be);
    }

    /// Solve the equilibrium problem
    auto solve(ChemicalState& state, double T, double P, VectorConstRef be) -> EquilibriumResult
    {
        // Check the dimension of the vector `be`
        Assert(be.size() == static_cast<int>(Ee),
            "Cannot proceed with method EquilibriumSolver::solve.",
            "The dimension of the given vector of molar amounts of the "
            "elements does not match the number of elements in the "
            "equilibrium partition.");
        return solve(state, T, P, be.data());
    }

    /// Solve the equilibrium problem
    auto solve(ChemicalState& state, double T, double P, const double* _be) -> EquilibriumResult
    {
        tic(SOLVE);

        // Reset the result of last equilibrium calculation
        result = {};

        // Set the molar amounts of the elements
        be = Vector::Map(_be, Ee);

        /*
        // If be is not changed from be_prev too much, skip this simulation
        Vector be_prev = state.elementAmountsInSpecies(ifs)(iee) + state.elementAmountsInSpecies(iss)(iee);
        const double diff = (be_prev - be).norm() / be.norm();

        if(diff < 1e-12)
        {
            skipped++;
            return result;
        }
        */
        // Set temperature and pressure of the chemical state
        state.setTemperature(T);
        state.setPressure(P);

        // Check if a simplex cold-start approximation must be performed
        if(coldstart(state))
            initialguess(state, T, P, be);

        // Update the optimum options
        updateOptimumOptions();

        // Update the optimum problem
        updateOptimumProblem(state);

        // Update the optimum state
        updateOptimumState(state);

        // Set the method for the optimisation calculation
        solver.setMethod(options.method);

        // Set the maximum number of iterations in each optimization pass
        optimum_options.max_iterations = 10;

        // Start the several opmization passes (stop if convergence attained)
        auto counter = 0;
        while(counter < options.optimum.max_iterations)
        {
            // Solve the optimisation problem
            result.optimum += solver.solve(optimum_problem, optimum_state, optimum_options);

            // Exit this loop if last solve succeeded
            if(result.optimum.succeeded)
                break;

            counter += optimum_options.max_iterations;
        }

        // Update the chemical state from the optimum state
        updateChemicalState(state);

        result.timing.solve = toc(SOLVE);

        return result;
    }

    /// Return the sensitivity of the equilibrium state.
    auto sensitivity() -> const EquilibriumSensitivity&
    {
        auto& dndT = sensitivities.dndT;
        auto& dydT = sensitivities.dydT;
        auto& dzdT = sensitivities.dzdT;
        auto& dndP = sensitivities.dndP;
        auto& dydP = sensitivities.dydP;
        auto& dzdP = sensitivities.dzdP;
        auto& dndb = sensitivities.dndb;
        auto& dydb = sensitivities.dydb;
        auto& dzdb = sensitivities.dzdb;

        auto& dnedT  = sensitivities_e.dndT;
        auto& dyedT  = sensitivities_e.dydT;
        auto& dzedT  = sensitivities_e.dzdT;
        auto& dnedP  = sensitivities_e.dndP;
        auto& dyedP  = sensitivities_e.dydP;
        auto& dzedP  = sensitivities_e.dzdP;
        auto& dnedbe = sensitivities_e.dndb;
        auto& dyedbe = sensitivities_e.dydb;
        auto& dzedbe = sensitivities_e.dzdb;

        const auto& dgedT  = ue.ddT;
        const auto& dgedP  = ue.ddP;
        const auto& dgedbe = zeros(Ne, Ee);
        const auto& dbedT  = zeros(Ee);
        const auto& dbedP  = zeros(Ee);
        const auto& dbedbe = identity(Ee, Ee);

        solver.sensitivities(dgedT,  dbedT,  dnedT,  dyedT,  dzedT);
        solver.sensitivities(dgedP,  dbedP,  dnedP,  dyedP,  dzedP);
        solver.sensitivities(dgedbe, dbedbe, dnedbe, dyedbe, dzedbe);

        dndT(ies)      = dnedT;
        dndP(ies)      = dnedP;
        dndb(ies, iee) = dnedbe;

        dydT(iee)      = dyedT;
        dydP(iee)      = dyedP;
        dydb(iee, iee) = dyedbe;

        dzdT(ies)      = dzedT;
        dzdP(ies)      = dzedP;
        dzdb(ies, iee) = dzedbe;

        return sensitivities;
    }
};

EquilibriumSolver::EquilibriumSolver()
{
    RuntimeError("Cannot proceed with EquilibriumSolver().",
        "EquilibriumSolver() constructor is deprecated. "
        "Use constructors EquilibriumSolver(const ChemicalSystem&) or EquilibriumSolver(const Partition&) instead.");
}

EquilibriumSolver::EquilibriumSolver(const ChemicalSystem& system)
: pimpl(new Impl(Partition(system)))
{}

EquilibriumSolver::EquilibriumSolver(const Partition& partition)
: pimpl(new Impl(partition))
{}

EquilibriumSolver::EquilibriumSolver(const EquilibriumSolver& other)
: pimpl(new Impl(*other.pimpl))
{}

EquilibriumSolver::~EquilibriumSolver()
{}

auto EquilibriumSolver::operator=(EquilibriumSolver other) -> EquilibriumSolver&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto EquilibriumSolver::setOptions(const EquilibriumOptions& options) -> void
{
    pimpl->options = options;
}

auto EquilibriumSolver::setPartition(const Partition& partition) -> void
{
    RuntimeError("Cannot proceed with EquilibriumSolver::setPartition.",
        "EquilibriumSolver::setPartition is deprecated. "
        "Use constructor EquilibriumSolver(const Partition&) instead.");
}

auto EquilibriumSolver::approximate(ChemicalState& state, double T, double P, VectorConstRef be) -> EquilibriumResult
{
    return pimpl->approximate(state, T, P, be);
}

auto EquilibriumSolver::approximate(ChemicalState& state, const EquilibriumProblem& problem) -> EquilibriumResult
{
    return pimpl->approximate_with_all_element_amounts(state, problem.temperature(), problem.pressure(), problem.elementAmounts());
}

auto EquilibriumSolver::approximate(ChemicalState& state) -> EquilibriumResult
{
    return pimpl->approximate_with_all_element_amounts(state, state.temperature(), state.pressure(), state.elementAmounts());
}

auto EquilibriumSolver::solve(ChemicalState& state) -> EquilibriumResult
{
    return pimpl->solve(state);
}

auto EquilibriumSolver::solve(ChemicalState& state, const EquilibriumProblem& problem) -> EquilibriumResult
{
    return pimpl->solve(state, problem);
}

auto EquilibriumSolver::solve(ChemicalState& state, double T, double P, VectorConstRef be) -> EquilibriumResult
{
    return pimpl->solve(state, T, P, be);
}

auto EquilibriumSolver::solve(ChemicalState& state, double T, double P, const double* be) -> EquilibriumResult
{
    return pimpl->solve(state, T, P, be);
}

auto EquilibriumSolver::properties() const -> const ChemicalProperties&
{
    return pimpl->properties;
}

auto EquilibriumSolver::sensitivity() -> const EquilibriumSensitivity&
{
    return pimpl->sensitivity();
}

auto EquilibriumSolver::result() const -> const EquilibriumResult&
{
    return pimpl->result;
}

} // namespace Reaktoro
