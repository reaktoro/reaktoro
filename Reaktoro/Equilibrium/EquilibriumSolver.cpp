// Reaktoro is a C++ library for computational reaction modelling.
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

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Math/MathUtils.hpp>
#include <Reaktoro/Optimization/OptimumOptions.hpp>
#include <Reaktoro/Optimization/OptimumProblem.hpp>
#include <Reaktoro/Optimization/OptimumResult.hpp>
#include <Reaktoro/Optimization/OptimumSolver.hpp>
#include <Reaktoro/Optimization/OptimumState.hpp>

namespace Reaktoro {
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

    /// The partition of the chemical system
    Partition partition;

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

    /// The indices of the equilibrium elements
    Indices iequilibrium_elements;

    /// The number of equilibrium species and elements
    unsigned Ne, Ee;

    /// The indices of the stable equilibrium species
    Indices istable_species;

    /// The indices of the unstable equilibrium species
    Indices iunstable_species;

    /// The formula matrix of the species
    Matrix A;

    /// The formula matrix of the stable species
    Matrix As;

    /// The formula matrix of the unstable species
    Matrix Au;

    /// Construct a Impl instance
    Impl(const ChemicalSystem& system)
    : system(system)
    {
        // Set the default partition as all species are in equilibrium
        setPartition(Partition(system));

        // Initialize the formula matrix
        A = system.formulaMatrix();
    }

    /// Set the partition of the chemical system
    auto setPartition(const Partition& partition_) -> void
    {
        // Set the partition of the chemical system
        partition = partition_;

        // Initialize the indices of the equilibrium species and elements
        iequilibrium_species = partition.indicesEquilibriumSpecies();
        iequilibrium_elements = partition.indicesEquilibriumElements();

        // Initialize the number of equilibrium species and elements
        Ne = iequilibrium_species.size();
        Ee = iequilibrium_elements.size();
    }

    /// Set the partition of the chemical system using a formatted string
    auto setPartition(std::string partition) -> void
    {
        setPartition(Partition(system, partition));
    }

    /// Update the sets of stable and unstable species by checking which species have zero molar amounts
    auto updateStabilitySets(const ChemicalState& state) -> void
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
        As = submatrix(A, iequilibrium_elements, istable_species);
        Au = submatrix(A, iequilibrium_elements, iunstable_species);
    }

    /// Update the molar amounts of unstable species based making them zero
    auto updateUnstableSpeciesAmounts(ChemicalState& state) -> void
    {
        for(Index i : istable_species)
            if(state.speciesAmount(i) < options.epsilon)
                state.setSpeciesAmount(i, 0.0);
    }

    /// Update the OptimumOptions instance with given EquilibriumOptions instance
    auto updateOptimumOptions() -> void
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
            for(Index i : iequilibrium_elements)
                ynames.push_back(system.element(i).name());

            // Initialize the names of the dual variables `z`
            znames = xnames;
        }
    }

    /// Update the OptimumProblem instance with given EquilibriumProblem and ChemicalState instances
    auto updateOptimumProblem(const ChemicalState& state, const Vector& be) -> void
    {
        // The number of stable species and elements in the equilibrium partition
        const unsigned num_stable_species = istable_species.size();
        const unsigned num_equilibrium_elements = iequilibrium_elements.size();

        // The temperature and pressure of the equilibrium calculation
        const double T  = state.temperature();
        const double P  = state.pressure();
        const double RT = universalGasConstant*T;

        // Set the molar amounts of the species
        n = state.speciesAmounts();

        // Set the Jacobian matrix of the optimisation calculation
        Jacobian jacobian;

        // The left-hand side matrix of the balance equations
        jacobian.Ae = As;

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
            return As*x - be;
        };

        // Define the gradient function of the mass-cahrge balance contraint function
        ConstraintGradFunction balance_constraint_grad = [=](const Vector& x) -> Jacobian
        {
            return jacobian;
        };

        // Setup an OptimumProblem instance with the Gibbs energy function and the balance constraints
        optimum_problem = OptimumProblem();
        optimum_problem.setNumVariables(num_stable_species);
        optimum_problem.setNumConstraints(num_equilibrium_elements);
        optimum_problem.setObjective(gibbs);
        optimum_problem.setObjectiveGrad(gibbs_grad);
        optimum_problem.setObjectiveHessian(gibbs_hessian);
        optimum_problem.setConstraint(balance_constraint);
        optimum_problem.setConstraintGrad(balance_constraint_grad);
        optimum_problem.setLowerBounds(0.0);
    }

    /// Initialize the optimum state from a chemical state
    auto updateOptimumState(const ChemicalState& state) -> void
    {
        // Set the molar amounts of the species
        n = state.speciesAmounts();

        // Set the dual potentials of the elements
        y = state.elementPotentials();

        // Set the dual potentials of the species
        z = state.speciesPotentials();

        // Initialize the optimum state
        rows(n, istable_species).to(optimum_state.x);
        rows(y, iequilibrium_elements).to(optimum_state.y);
        rows(z, istable_species).to(optimum_state.z);
    }

    /// Initialize the chemical state from a optimum state
    auto updateChemicalState(ChemicalState& state) -> void
    {
        // Update the dual potentials of the species and elements
        z.fill(0.0); rows(z, istable_species) = optimum_state.z;
        y.fill(0.0); rows(y, iequilibrium_elements) = optimum_state.y;

        // Update the chemical state
        state.setSpeciesAmounts(n);
        state.setElementPotentials(y);
        state.setSpeciesPotentials(z);
    }

    /// Find an initial feasible guess for an equilibrium problem
    auto approximate(ChemicalState& state, const Vector& be) -> EquilibriumResult
    {
        // Check the dimension of the vector `be`
        Assert(unsigned(be.rows()) == Ee,
            "Cannot proceed with method EquilibriumSolver::approximate.",
            "The dimension of the given vector of molar amounts of the "
            "elements does not match the number of elements in the "
            "equilibrium partition.");

        // The formula matrix of the equilibrium species
        const Matrix& Ae = partition.formulaMatrixEquilibriumSpecies();

        // The number of species and elements
        const unsigned num_species = system.numSpecies();
        const unsigned num_elements = system.numElements();

        // The number of equilibrium species and elements
        const unsigned num_equilibrium_species = iequilibrium_species.size();
        const unsigned num_equilibrium_elements = iequilibrium_elements.size();

        // The temperature and pressure of the equilibrium calculation
        const double T  = state.temperature();
        const double P  = state.pressure();
        const double RT = universalGasConstant*T;

        // Set the Jacobian matrix of the optimisation calculation
        Jacobian jacobian;

        // The left-hand side matrix of the mass balance equations
        jacobian.Ae = Ae;

        // Calculate the standard Gibbs energies of the species
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
            return Ae*ne - be;
        };

        // Define the gradient function of the mass-cahrge balance contraint function
        ConstraintGradFunction balance_constraint_grad = [=](const Vector& x) -> Jacobian
        {
            return jacobian;
        };

        // Setup an OptimumProblem instance with the Gibbs energy function and the balance constraints
        OptimumProblem optimum_problem;
        optimum_problem.setNumVariables(num_equilibrium_species);
        optimum_problem.setNumConstraints(num_equilibrium_elements);
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
        Vector y = zeros(num_elements);
        rows(z, iequilibrium_species) = optimum_state.z;
        rows(y, iequilibrium_elements)  = optimum_state.y;

        // Update the chemical state
        state.setSpeciesAmounts(n);
        state.setElementPotentials(y);
        state.setSpeciesPotentials(z);

        return result;
    }

    /// Solve the equilibrium problem
    auto solve(ChemicalState& state, const Vector& be) -> EquilibriumResult
    {
        // Check the dimension of the vector `be`
        Assert(unsigned(be.size()) == Ee,
            "Cannot proceed with method EquilibriumSolver::solve.",
            "The dimension of the given vector of molar amounts of the "
            "elements does not match the number of elements in the "
            "equilibrium partition.");

        // The result of the equilibrium calculation
        EquilibriumResult result;

        do {
            // Update the sets of stable and unstable species
            updateStabilitySets(state);

            // Update the optimum options
            updateOptimumOptions();

            // Update the optimum problem
            updateOptimumProblem(state, be);

            // Update the optimum state
            updateOptimumState(state);

            // Solve the optimisation problem
            result.optimum += solver.solve(optimum_problem, optimum_state, optimum_options);

            // Update the chemical state from the optimum state
            updateChemicalState(state);

        } while(not allStable(state));

        updateUnstableSpeciesAmounts(state);

        return result;
    }

    /// Return true if all species are stable
    auto allStable(ChemicalState& state) -> bool
    {
        if(iunstable_species.empty())
            return true;

        const double T = state.temperature();
        const double P = state.pressure();
        const double RT = universalGasConstant*T;
        const double lambda = std::sqrt(options.epsilon);

        for(Index i : iunstable_species)
            n[i] = options.epsilon;

        u = system.chemicalPotentials(T, P, n).val/RT;
        uu = rows(u, iunstable_species);
        zu = uu - tr(Au) * optimum_state.y;

        for(int k = 0; k < zu.rows(); ++k)
            if(zu[k] < 0.0) state.setSpeciesAmount(iunstable_species[k], lambda);

        return min(zu) >= 0.0;
    }

    /// Return the partial derivatives dn/db
    auto dndb(const ChemicalState& state) const -> Matrix
    {
        const unsigned Ne = iequilibrium_species.size();
        const unsigned Ee = iequilibrium_elements.size();

        const double& T = state.temperature();
        const double& P = state.pressure();
        const Vector& n = state.speciesAmounts();
        const Vector& y = state.elementPotentials();
        const Vector& z = state.speciesPotentials();

        const ChemicalVector u = system.chemicalPotentials(T, P, n);

        const Vector ne = rows(n, iequilibrium_species);
        const Vector ye = rows(y, iequilibrium_elements);
        const Vector ze = rows(z, iequilibrium_species);

        Jacobian Je;
        Je.Ae = submatrix(A, iequilibrium_elements, iequilibrium_species);

        Hessian He;
        He.mode = Hessian::Dense;
        He.dense = submatrix(u.ddn, iequilibrium_species, iequilibrium_species);

        KktSolution sol;

        KktMatrix lhs{He, Je, ne, ze};

        KktSolver kkt;
        kkt.decompose(lhs);

        KktVector rhs;
        rhs.rx = zeros(Ne);
        rhs.rz = zeros(Ne);

        Matrix dndb = zeros(Ne, Ee);

        for(Index i = 0; i < Ee; ++i)
        {
            rhs.ry = Vector::Unit(Ee, i);
            kkt.solve(rhs, sol);
            dndb.col(i) = sol.dx;
        }

        return dndb;
    }
};

EquilibriumSolver::EquilibriumSolver(const ChemicalSystem& system)
: pimpl(new Impl(system))
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
    pimpl->setPartition(partition);
}

auto EquilibriumSolver::setPartition(std::string partition) -> void
{
    pimpl->setPartition(partition);
}

auto EquilibriumSolver::approximate(ChemicalState& state, const Vector& be) -> EquilibriumResult
{
    return pimpl->approximate(state, be);
}

auto EquilibriumSolver::solve(ChemicalState& state, const Vector& be) -> EquilibriumResult
{
    // Solve the equilibrium problem
    EquilibriumResult result = pimpl->solve(state, be);

    // Return result if equilibrium calculation succeeded
    if(result.optimum.succeeded)
        return result;

    // Otherwise, solve it from scratch using an approximation
    result += pimpl->approximate(state, be);
    result += pimpl->solve(state, be);

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
    return pimpl->dndb(state);
}

} // namespace Reaktoro
