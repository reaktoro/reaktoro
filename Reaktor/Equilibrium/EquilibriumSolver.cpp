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
#include <Reaktor/Core/ChemicalState.hpp>
#include <Reaktor/Core/ChemicalSystem.hpp>
#include <Reaktor/Core/Partition.hpp>
#include <Reaktor/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktor/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktor/Equilibrium/EquilibriumResult.hpp>
#include <Reaktor/Optimization/OptimumOptions.hpp>
#include <Reaktor/Optimization/OptimumProblem.hpp>
#include <Reaktor/Optimization/OptimumResult.hpp>
#include <Reaktor/Optimization/OptimumSolver.hpp>
#include <Reaktor/Optimization/OptimumState.hpp>

namespace Reaktor {

struct EquilibriumSolver::Impl
{
    /// The options of the equilibrium solver
    EquilibriumOptions options;

    /// The solver for the optimisation calculations
    OptimumSolver solver;

    /// The molar amounts of the species
    Vector n;

    /// The dual potentials of the elements/charge
    Vector y;

    /// The dual potentials of the species
    Vector z;

    /// The molar amounts of the equilibrium species
    Vector ne;

    /// The chemical potentials of the species
    Vector u;

    /// The chemical potentials of the equilibrium species
    Vector ue;

    /// The optimisation problem
    OptimumProblem optimum_problem;

    /// The state of the optimisation calculation
    OptimumState optimum_state;

    /// Convert an EquilibriumProblem into an OptimumProblem
    auto updateOptimumProblem(const EquilibriumProblem& problem, const ChemicalState& state) -> void;

    /// Initialise the optimum state from a chemical state
    auto updateOptimumState(const EquilibriumProblem& problem, const ChemicalState& state) -> void;

    /// Initialise the chemical state from a optimum state
    auto updateChemicalState(const EquilibriumProblem& problem, ChemicalState& state) -> void;

    /// Find a initial guess for an equilibrium problem
    auto approximate(const EquilibriumProblem& problem, ChemicalState& state) -> EquilibriumResult;

    /// Solve the equilibrium problem
    auto solve(const EquilibriumProblem& problem, ChemicalState& state) -> EquilibriumResult;
};

auto EquilibriumSolver::Impl::updateOptimumProblem(const EquilibriumProblem& problem, const ChemicalState& state) -> void
{
    // The reference to the chemical system and its partitioning
    const ChemicalSystem& system = problem.system();
    const Partition& partition = problem.partition();

    // The indices of the equilibrium species
    const Indices& iequilibrium = partition.indicesEquilibriumSpecies();

    // The number of equilibrium species and linearly independent components in the equilibrium partition
    const unsigned num_equilibrium_species = iequilibrium.size();
    const unsigned num_components = problem.components().size();

    // The temperature and pressure of the equilibrium calculation
    const double T  = problem.temperature();
    const double P  = problem.pressure();
    const double RT = universalGasConstant*T;

    // Set the molar amounts of the species
    n = state.speciesAmounts();

    // Set the Jacobian matrix of the optimisation calculation
    Jacobian jacobian;

    // The left-hand side matrix of the linearly independent mass-charge balance equations
    jacobian.Ae = problem.balanceMatrix();

    // The right-hand side vector of the linearly independent mass-charge balance equations
    const Vector b = problem.componentAmounts();

    // Auxiliary function to update the state of a EquilibriumResult instance
    auto update = [=](const Vector& x) mutable
    {
        // Set the molar amounts of the species
        rows(n, iequilibrium) = x;

        // Set the molar amounts of the equilibrium species
        ne.noalias() = x;

        // Set the scaled chemical potentials of the species
        u = system.chemicalPotentials(T, P, n).val()/RT;

        // Set the scaled chemical potentials of the equilibrium species
        rows(u, iequilibrium).to(ue);
    };

    // Define the Gibbs energy function
    ObjectiveFunction gibbs = [=](const Vector& x) mutable
    {
        update(x);
        return dot(ne, ue);
    };

    // Define the gradient of the Gibbs energy function
    ObjectiveGradFunction gibbs_grad = [=](const Vector& x) mutable
    {
        if(x == ne) return ue;
        update(x);
        return ue;
    };

    // Define the Hessian of the Gibbs energy function
    ObjectiveHessianFunction gibbs_hessian = [=](const Vector& x, const Vector& g) mutable
    {
        Hessian hessian;
        hessian.mode = Hessian::Diagonal;
        hessian.diagonal = inv(x);

        // todo  todo  todo  todo  todo  todo  todo  todo  todo  todo  todo  todo
        for(unsigned i = 0; i < system.numPhases(); ++i)
            if(system.numSpeciesInPhase(i) == 1)
                hessian.diagonal[system.indexFirstSpeciesInPhase(i)] = 0.0;

        hessian.mode = Hessian::Dense;
        hessian.dense = diag(hessian.diagonal);

//        std::vector<double> nonzeros;
//        Indices inonzeros;
//        nonzeros.reserve(x.size());
//        inonzeros.reserve(x.size());
//
//        const auto& index_of_phase_with_species = system.connectivity().species_to_phase;
//
//        for(unsigned i = 0; i < system.numSpecies(); ++i)
//        {
//            const auto iphase = index_of_phase_with_species[i];
//            if(system.numSpeciesInPhase(iphase) > 1)
//            {
//                nonzeros.push_back(1.0/x[i]);
//                inonzeros.push_back(i);
//            }
//        }
//
//        hessian.mode = Hessian::SparseDiagonal;
//        hessian.sparsediagonal.nonzeros = Vector::Map(nonzeros.data(), nonzeros.size());
//        hessian.sparsediagonal.inonzeros = inonzeros;

        return hessian;
    };

    // Define the mass-cahrge balance contraint function
    ConstraintFunction balance_constraint = [=](const Vector& x) -> Vector
    {
        return jacobian.Ae*x - b;
    };

    // Define the gradient function of the mass-cahrge balance contraint function
    ConstraintGradFunction balance_constraint_grad = [=](const Vector& x) -> Jacobian
    {
        return jacobian;
    };

    // Setup an OptimumProblem instance with the Gibbs energy function and the balance constraints
    optimum_problem = OptimumProblem();
    optimum_problem.setNumVariables(num_equilibrium_species);
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
    // The reference to the chemical system
    const ChemicalSystem& system = problem.system();

    // The reference to the partition of the chemical system
    const Partition& partition = problem.partition();

    // The indices of the equilibrium species
    const Indices& iequilibrium = partition.indicesEquilibriumSpecies();

    // The indices of the linearly independent components (charge/elements)
    const Indices& icomponents = problem.components();

    // Set the molar amounts of the species
    n = state.speciesAmounts();

    // Set the dual potentials of the charge and elements
    y.conservativeResize(system.numElements() + 1);
    y << state.elementAmounts(), state.chargePotential();

    // Set the dual potentials of the species
    z = state.speciesPotentials();

    // Initialise the optimum state
    rows(n, iequilibrium).to(optimum_state.x);
    rows(y, icomponents).to(optimum_state.y);
    rows(z, iequilibrium).to(optimum_state.z);
}

auto EquilibriumSolver::Impl::updateChemicalState(const EquilibriumProblem& problem, ChemicalState& state) -> void
{
    // The reference to the chemical system
    const ChemicalSystem& system = problem.system();

    // The reference to the partition of the chemical system
    const Partition& partition = problem.partition();

    // The indices of the equilibrium species
    const Indices& iequilibrium = partition.indicesEquilibriumSpecies();

    // The indices of the linearly independent components (charge/elements)
    const Indices& icomponents = problem.components();

    // Update the dual potentials of the species and elements/charge
    z.fill(0.0); rows(z, iequilibrium) = optimum_state.z;
    y.fill(0.0); rows(y, icomponents)  = optimum_state.y;

    // Get the dual potential of electrical charge
    const double ycharge = y[system.numElements()];

    // Resize the dual potentials of charge/elements to remove the last entry corresponding to charge
    y.conservativeResize(system.numElements());

    // Update the chemical state
    state.setSpeciesAmounts(n);
    state.setChargePotential(ycharge);
    state.setElementPotentials(y);
    state.setSpeciesPotentials(z);
}

auto EquilibriumSolver::Impl::approximate(const EquilibriumProblem& problem, ChemicalState& state) -> EquilibriumResult
{
    // Initialise the optimum problem
    updateOptimumProblem(problem, state);

    // Initialise the optimum state
    updateOptimumState(problem, state);

    // The result of the equilibrium calculation
    EquilibriumResult result;

    // Find an approximation to the optimisation problem
    result.optimum = solver.approximate(optimum_problem, optimum_state, options.optimum);

    // Update the chemical state from the optimum state
    updateChemicalState(problem, state);

    return result;
}

auto EquilibriumSolver::Impl::solve(const EquilibriumProblem& problem, ChemicalState& state) -> EquilibriumResult
{
    // Initialise the optimum problem
    updateOptimumProblem(problem, state);

    // Initialise the optimum state
    updateOptimumState(problem, state);

    // The result of the equilibrium calculation
    EquilibriumResult result;

    // Solve the optimisation problem
    result.optimum = solver.solve(optimum_problem, optimum_state, options.optimum);

    // Update the chemical state from the optimum state
    updateChemicalState(problem, state);

    return result;
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
    return pimpl->solve(problem, state);
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
