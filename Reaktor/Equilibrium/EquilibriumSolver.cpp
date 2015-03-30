// Reaktor is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "EquilibriumSolver.hpp"

// Reaktor includes
#include <Reaktor/Common/Constants.hpp>
#include <Reaktor/Common/Exception.hpp>
#include <Reaktor/Core/ChemicalState.hpp>
#include <Reaktor/Core/ChemicalSystem.hpp>
#include <Reaktor/Core/Partition.hpp>
#include <Reaktor/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktor/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktor/Equilibrium/EquilibriumResult.hpp>
#include <Reaktor/Math/MathUtils.hpp>
#include <Reaktor/Optimization/OptimumOptions.hpp>
#include <Reaktor/Optimization/OptimumProblem.hpp>
#include <Reaktor/Optimization/OptimumResult.hpp>
#include <Reaktor/Optimization/OptimumSolver.hpp>
#include <Reaktor/Optimization/OptimumState.hpp>

namespace Reaktor {
namespace {

auto throwZeroInitialGuessError() -> void
{
    Exception exception;
    exception.error << "Cannot continue the equilibrium calculation.";
    exception.reason << "The provided initial state has zero molar amounts for all species.";
    RaiseError(exception);
}

} // namespace

struct EquilibriumSolver::Impl
{
    /// The chemical system instance
    ChemicalSystem system;

    /// The options of the equilibrium solver
    EquilibriumOptions options;

    /// The solver for the optimisation calculations
    OptimumSolver solver;

    /// The molar amounts of the species
    Vector n;

    /// The molar amounts of the stable equilibrium species
    Vector ns;

    /// The dual potentials of the elements
    Vector y;

    /// The dual potentials of the species
    Vector z;

    /// The dual potentials of the unstable species
    Vector zu;

    /// The chemical potentials of the species
    Vector u;

    /// The chemical potentials of the stable equilibrium species
    Vector us;

    /// The chemical potentials of the unstable equilibrium species
    Vector uu;

    /// The optimisation problem
    OptimumProblem optimum_problem;

    /// The state of the optimisation calculation
    OptimumState optimum_state;

    /// The options for the optimisation calculation
    OptimumOptions optimum_options;

    /// The indices of the equilibrium species
    Indices iequilibrium_species;

    /// The indices of the stable equilibrium species
    Indices istable_species;

    /// The indices of the unstable equilibrium species
    Indices iunstable_species;

    /// The indices of the stable elements
    Indices istable_elements;

    /// The formula matrix of the species
    Matrix A;

    /// The formula matrix of the stable species
    Matrix As;

    /// The formula matrix of the unstable species
    Matrix Au;

    /// Initialize the equilibrium solver
    auto initialize(const EquilibriumProblem& problem) -> void;

    /// Update the sets of stable and unstable species by checking which species have zero molar amounts
    auto updateStabilitySets(const ChemicalState& state) -> void;

    /// Update the molar amounts of unstable species based making them zero
    auto updateUnstableSpeciesAmounts(ChemicalState& state) -> void;

    /// Update the OptimumOptions instance with given EquilibriumOptions instance
    auto updateOptimumOptions() -> void;

    /// Update the OptimumProblem instance with given EquilibriumProblem and ChemicalState instances
    auto updateOptimumProblem(const EquilibriumProblem& problem, const ChemicalState& state) -> void;

    /// Initialize the optimum state from a chemical state
    auto updateOptimumState(const EquilibriumProblem& problem, const ChemicalState& state) -> void;

    /// Initialize the chemical state from a optimum state
    auto updateChemicalState(const EquilibriumProblem& problem, ChemicalState& state) -> void;

    /// Find an initial feasible guess for an equilibrium problem
    auto approximate(const EquilibriumProblem& problem, ChemicalState& state) -> EquilibriumResult;

    /// Solve the equilibrium problem
    auto solve(const EquilibriumProblem& problem, ChemicalState& state) -> EquilibriumResult;

    /// Return true if all species are stable
    auto allStable(const EquilibriumProblem& problem, ChemicalState& state) -> bool;
};

auto EquilibriumSolver::Impl::initialize(const EquilibriumProblem& problem) -> void
{
    // Initialize the chemical system
    system = problem.system();

    // Initialize the indices of the equilibrium species
    iequilibrium_species = problem.partition().indicesEquilibriumSpecies();

    // Initialize the formula matrix
    A = system.formulaMatrix();
}

auto EquilibriumSolver::Impl::updateStabilitySets(const ChemicalState& state) -> void
{
    // Update the indices of the stable and unstable species
    istable_species.clear();
    iunstable_species.clear();
    Vector n = state.speciesAmounts();
    for(Index i : iequilibrium_species)
        if(state.speciesAmount(i) > 0.0) istable_species.push_back(i);
        else iunstable_species.push_back(i);

    // Check if the set of stable species is empty
    if(istable_species.empty())
        throwZeroInitialGuessError();

    // Update the balance matrices of stable and unstable species
    As = cols(A, istable_species);
    Au = cols(A, iunstable_species);

    // Update the indices of the stable components
    istable_elements = linearlyIndependentRows(As, As);
}

auto EquilibriumSolver::Impl::updateUnstableSpeciesAmounts(ChemicalState& state) -> void
{
    for(Index i : istable_species)
        if(state.speciesAmount(i) < options.epsilon)
            state.setSpeciesAmount(i, 0.0);
}

auto EquilibriumSolver::Impl::updateOptimumOptions() -> void
{
    // Initialize the options for the optimisation calculation
    optimum_options = options.optimum;
    optimum_options.ipnewton.mu = options.epsilon * 1e-5;
    optimum_options.ipopt.mu.push_back(options.epsilon * 1e-5);

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
        for(Index i : istable_species)
            xnames.push_back(system.species(i).name());

        // Initialize the names of the dual variables `y`
        for(Index i : istable_elements)
            ynames.push_back(i < system.numElements() ? system.element(i).name() : "Z");

        // Initialize the names of the dual variables `z`
        znames = xnames;
    }
}

auto EquilibriumSolver::Impl::updateOptimumProblem(const EquilibriumProblem& problem, const ChemicalState& state) -> void
{
    // The number of stable equilibrium species and linearly independent components in the equilibrium partition
    const unsigned num_stable_species = istable_species.size();
    const unsigned num_components = istable_elements.size();

    // The temperature and pressure of the equilibrium calculation
    const double T  = problem.temperature();
    const double P  = problem.pressure();
    const double RT = universalGasConstant*T;

    // Set the molar amounts of the species
    n = state.speciesAmounts();

    // Set the Jacobian matrix of the optimisation calculation
    Jacobian jacobian;

    // The left-hand side matrix of the balance equations
    jacobian.Ae = As;

    // The right-hand side vector of the balance equations
    const Vector b = problem.elementAmounts();

    // Extract the stable components only
    Vector bs = rows(b, istable_elements);
    bs = max(bs, options.epsilon*ones(bs.rows()) * 1e-5);

    // Auxiliary function to update the state of a EquilibriumResult instance
    auto update = [=](const Vector& x) mutable
    {
        // Set the molar amounts of the species
        rows(n, istable_species) = x;

        // Set the molar amounts of the equilibrium species
        ns.noalias() = x;

        // Set the scaled chemical potentials of the species
        u = system.chemicalPotentials(T, P, n).val/RT;

        // Set the scaled chemical potentials of the equilibrium species
        rows(u, istable_species).to(us);
    };

    // Define the Gibbs energy function
    ObjectiveFunction gibbs = [=](const Vector& x) mutable
    {
        update(x);
        return dot(ns, us);
    };

    // Define the gradient of the Gibbs energy function
    ObjectiveGradFunction gibbs_grad = [=](const Vector& x) mutable
    {
        if(x == ns) return us;
        update(x);
        return us;
    };

    // Define the Hessian of the Gibbs energy function
    ObjectiveHessianFunction gibbs_hessian = [=](const Vector& x, const Vector& g) mutable
    {
        Hessian hessian;
        hessian.mode = Hessian::Diagonal;
        hessian.diagonal = inv(x);

        for(Index i = 0; i < istable_species.size(); ++i)
            if(system.numSpeciesInPhase(
                system.indexPhaseWithSpecies(istable_species[i])) == 1)
                    hessian.diagonal[i] = 0.0;

        return hessian;
    };

    // Define the mass-cahrge balance contraint function
    ConstraintFunction balance_constraint = [=](const Vector& x) -> Vector
    {
        return As*x - bs;
    };

    // Define the gradient function of the mass-cahrge balance contraint function
    ConstraintGradFunction balance_constraint_grad = [=](const Vector& x) -> Jacobian
    {
        return jacobian;
    };

    // Setup an OptimumProblem instance with the Gibbs energy function and the balance constraints
    optimum_problem = OptimumProblem();
    optimum_problem.setNumVariables(num_stable_species);
    optimum_problem.setNumConstraints(num_components);
    optimum_problem.setObjective(gibbs);
    optimum_problem.setObjectiveGrad(gibbs_grad);
    optimum_problem.setObjectiveHessian(gibbs_hessian);
    optimum_problem.setConstraint(balance_constraint);
    optimum_problem.setConstraintGrad(balance_constraint_grad);
    optimum_problem.setLowerBounds(0.0);
}

auto EquilibriumSolver::Impl::updateOptimumState(const EquilibriumProblem& problem, const ChemicalState& state) -> void
{
    // Set the molar amounts of the species
    n = state.speciesAmounts();

    // Set the dual potentials of the elements
    y = state.elementPotentials();

    // Set the dual potentials of the species
    z = state.speciesPotentials();

    // Initialize the optimum state
    rows(n, istable_species).to(optimum_state.x);
    rows(y, istable_elements).to(optimum_state.y);
    rows(z, istable_species).to(optimum_state.z);
}

auto EquilibriumSolver::Impl::updateChemicalState(const EquilibriumProblem& problem, ChemicalState& state) -> void
{
    // Update the dual potentials of the species and elements
    z.fill(0.0); rows(z, istable_species) = optimum_state.z;
    y.fill(0.0); rows(y, istable_elements)  = optimum_state.y;

    // Update the chemical state
    state.setSpeciesAmounts(n);
    state.setElementPotentials(y);
    state.setSpeciesPotentials(z);
}

auto EquilibriumSolver::Impl::approximate(const EquilibriumProblem& problem, ChemicalState& state) -> EquilibriumResult
{
    // Get a reference to the chemical system and its partitioning
    const ChemicalSystem& system = problem.system();
    const Partition& partition = problem.partition();

    // The indices of the equilibrium species
    const Indices& iequilibrium_species = partition.indicesEquilibriumSpecies();

    // Assemble the balance matrix of the equilibrium species only
    const Matrix Ae = cols(A, iequilibrium_species);

    // The indices of the linearly independent elements
    const Indices icomponents = linearlyIndependentRows(Ae);

    // Assemble the balance matrix of with stable elements
    const Matrix As = rows(Ae, icomponents);

    // The number of equilibrium species and linearly independent components in the equilibrium partition
    const unsigned num_species = system.numSpecies();
    const unsigned num_equilibrium_species = iequilibrium_species.size();
    const unsigned num_elements = system.numElements();
    const unsigned num_components = icomponents.size();

    // The temperature and pressure of the equilibrium calculation
    const double T  = problem.temperature();
    const double P  = problem.pressure();
    const double RT = universalGasConstant*T;

    // Set the Jacobian matrix of the optimisation calculation
    Jacobian jacobian;

    // The left-hand side matrix of the linearly independent balance equations
    jacobian.Ae = As;

    // The right-hand side vector of the linearly independent balance equations
    const Vector b = problem.elementAmounts();
    const Vector bs = rows(b, icomponents);

    const Vector u0 = system.standardGibbsEnergies(T, P).val/RT;
    const Vector ue0 = rows(u0, iequilibrium_species);

    Hessian hessian;
    hessian.mode = Hessian::Diagonal;
    hessian.diagonal = zeros(num_equilibrium_species);

    // Define the Gibbs energy function
    ObjectiveFunction gibbs = [=](const Vector& ne) mutable
    {
        return dot(ne, ue0);
    };

    // Define the gradient of the Gibbs energy function
    ObjectiveGradFunction gibbs_grad = [=](const Vector& ne) mutable
    {
        return ue0;
    };

    // Define the Hessian of the Gibbs energy function
    ObjectiveHessianFunction gibbs_hessian = [=](const Vector& x, const Vector& g) mutable
    {
        return hessian;
    };

    // Define the mass-cahrge balance contraint function
    ConstraintFunction balance_constraint = [=](const Vector& ne) -> Vector
    {
        return As*ne - bs;
    };

    // Define the gradient function of the mass-cahrge balance contraint function
    ConstraintGradFunction balance_constraint_grad = [=](const Vector& x) -> Jacobian
    {
        return jacobian;
    };

    // Setup an OptimumProblem instance with the Gibbs energy function and the balance constraints
    OptimumProblem optimum_problem;
    optimum_problem.setNumVariables(num_equilibrium_species);
    optimum_problem.setNumConstraints(num_components);
    optimum_problem.setObjective(gibbs);
    optimum_problem.setObjectiveGrad(gibbs_grad);
    optimum_problem.setObjectiveHessian(gibbs_hessian);
    optimum_problem.setConstraint(balance_constraint);
    optimum_problem.setConstraintGrad(balance_constraint_grad);
    optimum_problem.setLowerBounds(0.0);

    Vector n = state.speciesAmounts();

    // Initialize the options for the optimisation calculation
    OptimumOptions optimum_options = options.optimum;
    optimum_options.ipnewton.mu = options.epsilon * 1e-5;
    optimum_options.ipopt.mu.push_back(options.epsilon * 1e-5);

    // Initialize the optimum state
    OptimumState optimum_state;
    rows(n, iequilibrium_species).to(optimum_state.x);

    // The result of the equilibrium calculation
    EquilibriumResult result;

    // Find an approximation to the optimisation problem
    result.optimum = solver.solve(optimum_problem, optimum_state, optimum_options);

    // Update the chemical state from the optimum state
    rows(n, iequilibrium_species) = optimum_state.x;

    // Update the dual potentials of the species and elements
    Vector z = zeros(num_species);
    Vector y = zeros(num_elements + 1);
    rows(z, iequilibrium_species) = optimum_state.z;
    rows(y, icomponents)  = optimum_state.y;

    // Update the chemical state
    state.setSpeciesAmounts(n);
    state.setElementPotentials(y);
    state.setSpeciesPotentials(z);

    return result;
}

auto EquilibriumSolver::Impl::solve(const EquilibriumProblem& problem, ChemicalState& state) -> EquilibriumResult
{
    // The result of the equilibrium calculation
    EquilibriumResult result;

    // Initialize the equilibrium solver
    initialize(problem);

    do {
        // Update the sets of stable and unstable species
        updateStabilitySets(state);

        // Update the optimum options
        updateOptimumOptions();

        // Update the optimum problem
        updateOptimumProblem(problem, state);

        // Update the optimum state
        updateOptimumState(problem, state);

        // Solve the optimisation problem
        result.optimum += solver.solve(optimum_problem, optimum_state, optimum_options);

        // Update the chemical state from the optimum state
        updateChemicalState(problem, state);

    } while(not allStable(problem, state));

    updateUnstableSpeciesAmounts(state);

    return result;
}

auto EquilibriumSolver::Impl::allStable(const EquilibriumProblem& problem, ChemicalState& state) -> bool
{
    if(iunstable_species.empty())
        return true;

    const double T = problem.temperature();
    const double P = problem.pressure();
    const double RT = universalGasConstant*T;
    const double lambda = std::sqrt(options.epsilon);

    for(Index i : iunstable_species)
        n[i] = options.epsilon;

    u = system.chemicalPotentials(T, P, n).val/RT;
    uu = rows(u, iunstable_species);
    zu = uu - tr(Au) * y;

    for(int k = 0; k < zu.rows(); ++k)
        if(zu[k] < 0.0) state.setSpeciesAmount(iunstable_species[k], lambda);

    return min(zu) >= 0.0;
}

EquilibriumSolver::EquilibriumSolver()
: pimpl(new Impl())
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

auto EquilibriumSolver::approximate(const EquilibriumProblem& problem, ChemicalState& state) -> EquilibriumResult
{
    return pimpl->approximate(problem, state);
}

auto EquilibriumSolver::solve(const EquilibriumProblem& problem, ChemicalState& state) -> EquilibriumResult
{
    // Update the temperature and pressure of the chemical state
    state.setTemperature(problem.temperature());
    state.setPressure(problem.pressure());

    // Solve the equilibrium problem
    EquilibriumResult result = pimpl->solve(problem, state);

    // Return result if equilibrium calculation succeeded
    if(result.optimum.succeeded)
        return result;

    // Otherwise, solve it from scratch
    state.setSpeciesAmounts(0.0);
    result += pimpl->approximate(problem, state);
    result += pimpl->solve(problem, state);

    return result;
}

auto EquilibriumSolver::dndt(const ChemicalState& state) -> Vector
{
    return Vector();
}

auto EquilibriumSolver::dndp(const ChemicalState& state) -> Vector
{
    return Vector();
}

auto EquilibriumSolver::dndb(const ChemicalState& state) -> Matrix
{
    return Matrix();
}

} // namespace Reaktor
