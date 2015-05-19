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
    ChemicalVector u;

    /// The chemical potentials of the stable equilibrium species
    ChemicalVector us;

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

    /// The flags that indicate if an element always exist in positive amounts
    std::vector<bool> is_positive_element;

    /// The number of equilibrium species and elements
    unsigned Ne, Ee;

    /// The indices of the stable equilibrium species
    Indices istable_species;

    /// The indices of the unstable equilibrium species
    Indices iunstable_species;

    /// The formula matrix of the species
    Matrix A;

    /// The formula matrix of the equilibrium species
    Matrix Ae;

    /// The formula matrix of the stable species
    Matrix As;

    /// The formula matrix of the unstable species
    Matrix Au;

    /// Construct a Impl instance
    Impl(const ChemicalSystem& system)
    : system(system)
    {
        // Initialize the formula matrix
        A = system.formulaMatrix();

        // Set the default partition as all species are in equilibrium
        setPartition(Partition(system));

        // Determine the indices of the elements that always exist in positive amounts
        is_positive_element.resize(A.rows(), false);
        for(int i = 0; i < A.rows(); ++i)
            if(min(A.row(i)) >= 0)
                is_positive_element[i] = true;
    }

    /// Set the partition of the chemical system
    auto setPartition(const Partition& partition_) -> void
    {
        // Set the partition of the chemical system
        partition = partition_;

        // Initialize the formula matrix of the equilibrium species
        Ae = partition.formulaMatrixEquilibriumSpecies();

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

    // Enforce positive molar amounts for positive elements.
    auto regularizeElementAmounts(Vector& be) -> void
    {
        // Check if each equilibrium element has positive amounts
        for(unsigned i = 0; i < Ee; ++i)
        {
            const auto index = iequilibrium_elements[i];
            if(is_positive_element[index] and be[i] <= 0)
                be[i] = options.epsilon;
        }
    }

    /// Update the sets of stable and unstable species by checking which species have zero molar amounts
    auto updateStabilitySets(const ChemicalState& state) -> void
    {
        // Update the indices of the stable and unstable species
        istable_species.clear();
        iunstable_species.clear();
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

        // The temperature and pressure of the equilibrium calculation
        const double T  = state.temperature();
        const double P  = state.pressure();
        const double RT = universalGasConstant*T;

        // Set the molar amounts of the species
        n = state.speciesAmounts();

        // The result of the objective evaluation
        ObjectiveResult res;

        optimum_problem.objective = [=](const Vector& x) mutable
        {
            // Set the molar amounts of the species
            rows(n, istable_species) = x;

            // Set the molar amounts of the equilibrium species
            ns.noalias() = x;

            // Set the scaled chemical potentials of the species
            u = system.chemicalPotentials(T, P, n)/RT;

            // Set the scaled chemical potentials of the equilibrium species
            us = u.rows(istable_species, istable_species);

            // Set the objective result
            res.val = dot(ns, us.val);
            res.grad = us.val;

            if(options.hessian == EquilibriumHessian::Diagonal)
            {
                res.hessian.mode = Hessian::Diagonal;
                res.hessian.diagonal = inv(x);

                for(Index i = 0; i < istable_species.size(); ++i)
                    if(system.numSpeciesInPhase(
                        system.indexPhaseWithSpecies(istable_species[i])) == 1)
                            res.hessian.diagonal[i] = 0.0;
            }
            else
            {
                res.hessian.mode = Hessian::Dense;
                res.hessian.dense = us.ddn;
            }

            return res;
        };

        optimum_problem.A = As;
        optimum_problem.b = be;
        optimum_problem.l = zeros(num_stable_species);
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
    auto approximate(ChemicalState& state, Vector be) -> EquilibriumResult
    {
        // Check the dimension of the vector `be`
        Assert(unsigned(be.rows()) == Ee,
            "Cannot proceed with method EquilibriumSolver::approximate.",
            "The dimension of the given vector of molar amounts of the "
            "elements does not match the number of elements in the "
            "equilibrium partition.");

        // Enforce positive molar amounts for positive elements
        regularizeElementAmounts(be);

        // The temperature and pressure of the equilibrium calculation
        const double T  = state.temperature();
        const double P  = state.pressure();
        const double RT = universalGasConstant*T;

        // Calculate the standard Gibbs energies of the species
        const Vector u0 = system.standardGibbsEnergies(T, P).val/RT;
        const Vector ue0 = rows(u0, iequilibrium_species);

        ObjectiveResult res;
        res.hessian.mode = Hessian::Diagonal;
        res.hessian.diagonal = zeros(Ne);

        // Define the optimisation problem
        OptimumProblem optimum_problem;
        optimum_problem.A = Ae;
        optimum_problem.b = be;
        optimum_problem.l = zeros(Ne);
        optimum_problem.objective = [=](const Vector& ne) mutable
        {
            res.val = dot(ne, ue0);
            res.grad = ue0;
            return res;
        };

        n = state.speciesAmounts();
        y = state.elementPotentials();
        z = state.speciesPotentials();

        // Initialize the options for the optimisation calculation
        OptimumOptions optimum_options = options.optimum;
        optimum_options.kkt.method = KktMethod::Rangespace;
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
            for(Index i : iequilibrium_species)
                xnames.push_back(system.species(i).name());

            // Initialize the names of the dual variables `y`
            for(Index i : iequilibrium_elements)
                ynames.push_back(system.element(i).name());

            // Initialize the names of the dual variables `z`
            znames = xnames;
        }

        // Initialize the optimum state
        OptimumState optimum_state;
        rows(n, iequilibrium_species).to(optimum_state.x);
        rows(y, iequilibrium_elements).to(optimum_state.y);
        rows(z, iequilibrium_species).to(optimum_state.z);

        // The result of the equilibrium calculation
        EquilibriumResult result;

        // Find an approximation to the optimisation problem
        result.optimum = solver.solve(optimum_problem, optimum_state, optimum_options);

        // Update the chemical state from the optimum state
        rows(n, iequilibrium_species) = optimum_state.x;

        // Update the dual potentials of the species and elements
        z.fill(0.0); rows(z, iequilibrium_species) = optimum_state.z;
        y.fill(0.0); rows(y, iequilibrium_elements)  = optimum_state.y;

        // Update the chemical state
        state.setSpeciesAmounts(n);
        state.setElementPotentials(y);
        state.setSpeciesPotentials(z);

        return result;
    }

    /// Solve the equilibrium problem
    auto solve(ChemicalState& state, Vector be) -> EquilibriumResult
    {
        // Check the dimension of the vector `be`
        Assert(unsigned(be.size()) == Ee,
            "Cannot proceed with method EquilibriumSolver::solve.",
            "The dimension of the given vector of molar amounts of the "
            "elements does not match the number of elements in the "
            "equilibrium partition.");

        // Enforce positive molar amounts for positive elements
        regularizeElementAmounts(be);

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

            // Exit if the the optimisation calculation did not succeed
            if(not result.optimum.succeeded) break;

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

        u = system.chemicalPotentials(T, P, n)/RT;
        uu = rows(u.val, iunstable_species);
        zu = uu - tr(Au) * optimum_state.y;

        for(int k = 0; k < zu.rows(); ++k)
            if(zu[k] < 0.0) state.setSpeciesAmount(iunstable_species[k], lambda);

        return min(zu) >= 0.0;
    }

    /// Return the partial derivatives dn/db
    auto dndb(const ChemicalState& state) -> Matrix
    {
        updateStabilitySets(state);

        const Vector ns = rows(n, istable_species);
        const Vector zs = rows(z, istable_species);

        Hessian Hs;
        Hs.mode = Hessian::Dense;
        Hs.dense = submatrix(u.ddn, istable_species, istable_species);

        KktSolution sol;

        Matrix Ass;
        const Indices ielements = linearlyIndependentRows(As, Ass);

        KktMatrix lhs{Hs, Ass, ns, zs};

        KktSolver kkt;
        kkt.decompose(lhs);

        const unsigned Ns = istable_species.size();

        KktVector rhs;
        rhs.rx = zeros(Ns);
        rhs.rz = zeros(Ns);

        const unsigned Es = ielements.size();

        Matrix dnsdbs = zeros(Ns, Es);

        for(Index i = 0; i < ielements.size(); ++i)
        {
            rhs.ry = Vector::Unit(Es, i);
            kkt.solve(rhs, sol);
            dnsdbs.col(i) = sol.dx;
        }

        Matrix dndb = zeros(Ne, Ee);
        submatrix(dndb, istable_species, ielements) = dnsdbs;

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
