/*
 * Reaktor is a C++ library for computational reaction modelling.
 *
 * Copyright (C) 2014 Allan Leal
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "EquilibriumSolver.hpp"

// C++ includes
#include <cmath>

// Eigen includes
#include <Eigen/Dense>

// Reaktor includes
#include <Reaktor/Common/Constants.hpp>
#include <Reaktor/Common/Exception.hpp>
#include <Reaktor/Core/ChemicalState.hpp>
#include <Reaktor/Core/ChemicalSystem.hpp>
#include <Reaktor/Core/Partitioning.hpp>
#include <Reaktor/Equilibrium/BalanceConstraints.hpp>
#include <Reaktor/Equilibrium/EquilibriumConstraints.hpp>
#include <Reaktor/Equilibrium/EquilibriumLagrange.hpp>
#include <Reaktor/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktor/Equilibrium/EquilibriumParams.hpp>
#include <Reaktor/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktor/Equilibrium/EquilibriumResult.hpp>
#include <Reaktor/Utils/MatrixUtils.hpp>

namespace Reaktor {
namespace internal {

auto equilibriumObjective(const ChemicalSystem& system, const Partitioning& partitioning, ChemicalState state) -> Optima::ObjectiveFunction
{
    // The temperature of the chemical system
    const double T = state.temperature();

    // The pressure of the chemical system
    const double P = state.pressure();

    // The universal gas constant times temperature variable
    const double RT = universalGasConstant * T;

    // The standard chemical potentials of the equilibrium species
    const Vector ue0 = partitioning.equilibriumRows(system.chemicalPotentials(T, P));

    // The activities of the chemical species
    VectorResult a;

    // The activities of the equilibrium species
    VectorResult ae;

    // The natural logarithm of the activities of the equilibrium species
    Vector ln_ae;

    // The chemical potentials of the equilibrium species
    Vector ue;

    // The partial molar derivatives of the natural logarithm of the activities of the equilibrium species
    Matrix Ae;

    // The result of the objective function evaluation
    Optima::ObjectiveResult result;

    // The objective function
    Optima::ObjectiveFunction fn = [=,&partitioning](const Vector& ne) mutable
    {
        // Set the components of `n` that corresponds to the equilibrium species
        state.setEquilibriumComposition(partitioning, ne);

        // Calculate the activities of the species and their derivatives
        a = state.activities();

        // Extract the rows of `a` that corresponds to the equilibrium species
        func(ae) = partitioning.equilibriumRows(func(a));

        // Extract the rows and cols of `A` that corresponds to the equilibrium species
        grad(ae) = partitioning.equilibriumRowsCols(grad(a));

        // Calculate the partial molar derivatives of the natural logarithm of the activities of the equilibrium species
        Ae = func(ae).asDiagonal().inverse() * grad(ae);

        // Calculate the natural logarithm of the activities of the equilibrium species
        ln_ae = func(ae).array().log();

        // Calculate the normalised chemical potentials of the equilibrium species
        ue = ue0/RT + ln_ae;

        // Calculate the Gibbs free energy of the equilibrium system
        result.func = ue.dot(ne);

        // Calculate the first order partial derivatives of Gibbs free energy of the equilibrium system
        result.grad = ue;

        // Calculate the second order partial derivatives of Gibbs free energy of the equilibrium system
        result.hessian = Ae;

        return result;
    };

    return fn;
}

Optima::ConstraintFunction equilibriumConstraint(const EquilibriumConstraints& constraint, const Partitioning& partitioning, ChemicalState state)
{
    // The result of a vector function evaluation
    VectorResult he;

    // The result of the constraint function evaluation
    Optima::ConstraintResult result;

    // The constraint function
    Optima::ConstraintFunction fn = [=](const Vector& ne) mutable
    {
        state.setEquilibriumComposition(partitioning, ne);
        he = constraint(state);
        result.func = func(he);
        result.grad = grad(he);

        return result;
    };

    return fn;
}

EquilibriumConstraints massBalance(const ChemicalSystem& system, const Partitioning& partitioning, const Vector& be)
{
    Matrix We = partitioning.equilibriumFormulaMatrix(system);
    Vector ne;
    VectorResult he;
    grad(he) = We;

    auto fn = [=](const ChemicalState& state) mutable
    {
        ne = state.equilibriumComposition(partitioning);
        func(he) = We * ne - be;
        return he;
    };

    return EquilibriumConstraints(fn, We.rows());
}

} /* namespace internal */

double calculateLargestBoundaryStep(const Vector& p, const Vector& dp)
{
    double step = 1.0e+300;
    for(unsigned i = 0; i < p.rows(); ++i)
    {
        const double aux = -p[i]/dp[i];
        if(aux > 0.0) step = std::min(step, aux);
    }
    return (step > 1.0) ? 1.0 : 0.995 * step;
}

using namespace internal;

class EquilibriumSolver::Impl
{
    /// The chemical system instance
    ChemicalSystem system$;

    /// The partitioning of the species in the chemical system
    Partitioning partitioning$;

    /// The mass-balance and charge-balance constraints
    BalanceConstraints balance$;

    /// The equilibrium options of the calculation
    EquilibriumOptions options$;

    /// The equilibrium params of the calculation
    EquilibriumParams params$;

    /// The minimisation solver used to minimise the Gibbs free energy of the system
    Optima::IPFilterSolver minimum_solver$;

    /// The scaling to be performed in the equilibrium calculation
    Optima::Scaling scaling$;

    /// The formula matrix of the equilibrium species
    Matrix We$;

public:
    Impl(const ChemicalSystem& system, const Partitioning& partitioning)
    : system$(system), partitioning$(partitioning),
      balance$(system, partitioning),
      We$(partitioning$.equilibriumFormulaMatrix(system$))
    {
        minimum_solver$.SetParams(EquilibriumParams());
    }

    auto setOptions(const EquilibriumOptions& options) -> void
    {
        options$ = options;
        minimum_solver$.SetOptions(options);
    }

    auto setParams(const EquilibriumParams& params) -> void
    {
        params$ = params;
        minimum_solver$.SetParams(params);
    }

    auto setScaling(const Optima::Scaling& scaling) -> void
    {
        scaling$ = scaling;
    }

    auto setScaling(const ChemicalState& state) -> void
    {
        // Extract the molar abundance of the equilibrium species
        const Vector ne = state.equilibriumComposition(partitioning$);

        // Set the scaling of the variables
        scaling$.SetScalingVariables(ne);
    }

    auto minimise(ChemicalState& state, EquilibriumLagrange& lagrange, const EquilibriumConstraints& constraints) -> EquilibriumResult
    {
        // Update the scaling in the minimisation solver
        minimum_solver$.SetScaling(scaling$);

        // The number of elements and equilibrium species
        const unsigned Nc = constraints.numConstraints();
        const unsigned Ne = partitioning$.numEquilibriumSpecies();

        // Creates the objective function instance
        const auto objectivefn = equilibriumObjective(system$, partitioning$, state);

        // Creates the constraint function instance
        const auto constraintfn = equilibriumConstraint(constraints, partitioning$, state);

        // Creates the optimisation problem
        Optima::OptimumProblem opt_problem;
        opt_problem.SetNumVariables(Ne);
        opt_problem.SetNumConstraints(Nc);
        opt_problem.SetObjectiveFunction(objectivefn);
        opt_problem.SetConstraintFunction(constraintfn);

        // Set the optimisation problem in the minimisation solver
        minimum_solver$.SetProblem(opt_problem);

        // Extract the molar composition of the equilibrium species
        Vector ne = state.equilibriumComposition(partitioning$);

        // The normalised Lagrange multipliers of the equilibrium calculation
        Vector ye = lagrange.lagrangeY();
        Vector ze = lagrange.lagrangeZ();

        // Solve the minimisation problem
        minimum_solver$.Solve(ne, ye, ze);

        // Set the molar composition of the equilibrium species
        state.setEquilibriumComposition(partitioning$, ne);

        // Set the Lagrange multipliers of the equilibrium problem
        lagrange.setMultipliers(state, ye, ze);

        return minimum_solver$.GetResult();
    }

    auto minimise(ChemicalState& state, EquilibriumLagrange& lagrange, const Vector& be) -> EquilibriumResult
    {
        balance$.setMassBalance(be);
        return minimise(state, lagrange, balance$);
    }

    auto refine(ChemicalState& state, EquilibriumLagrange& lagrange) -> EquilibriumResult
    {
        EquilibriumResult result;

        // The molar abundance of the equilibrium species
        const Vector ne = state.equilibriumComposition(partitioning$);

        // The molar abundance of the equilibrium elements
        const Vector be = We$ * ne;

        // Set the molar abundance of the equilibrium elements in the mass-charge-balance constraint
        balance$.setMassBalance(be);

        // The temperature of the chemical system
        const double T = state.temperature();

        // The pressure of the chemical system
        const double P = state.pressure();

        // The temperature times universal gas constant parameter
        const double RT = universalGasConstant * T;

        // The indices of the stable equilibrium species (indices relative to the set of all species)
        const Indices& idx_stable_species = lagrange.idxStableSpecies();

        // The indices of the unstable equilibrium species (indices relative to the set of all species)
        const Indices& idx_unstable_species = lagrange.idxUnstableSpecies();

        // The indices of the stable equilibrium species (indices relative to the set of equilibrium species)
        const Indices& idx_equilibrium_stable_species = lagrange.idxEquilibriumStableSpecies();

        // The indices of the unstable equilibrium species (indices relative to the set of equilibrium species)
        const Indices& idx_equilibrium_unstable_species = lagrange.idxEquilibriumUnstableSpecies();

        // The mass and charge balance matrix of the equilibrium species
        const Matrix& He = balance$.balanceMatrix();

        // Initialise the balance matrix of the stable equilibrium species
        const Matrix Hs = cols(idx_equilibrium_stable_species, He);

        // Initialise the balance matrix of the unstable equilibrium species
        const Matrix Hu = cols(idx_equilibrium_unstable_species, He);

        // The number of equilibrium mass-charge balance constraints
        const unsigned Nc = balance$.numConstraints();

        // The number of stable equilibrium species
        const unsigned Ns = idx_stable_species.size();

        // The number of equilibrium reactions among the stable species
        const unsigned Ms = Ns - Nc;

        // Initialise the kernel of the formula matrix of the stable species
        Matrix vs = Hs.fullPivLu().kernel().transpose();

        // Calculate the standard chemical potential of the species
        const Vector u0 = system$.chemicalPotentials(T, P);

        // Calculate the standard chemical potential of the stable species
        const Vector us0 = rows(idx_stable_species, u0);

        // Calculate the standard chemical potential of the unstable species
        const Vector uu0 = rows(idx_unstable_species, u0);

        // Calculate the equilibrium constants of the reactions among the stable species
        const Vector ln_Ks = -1.0/RT * vs * us0;

        // Set the molar composition of the unstable species to zero
        state.setComposition(idx_unstable_species, Vector::Constant(Ns, 1.0e-50));

        // The molar composition of the stable equilibrium species
        Vector ns = rows(idx_stable_species, state.composition());

        // Auxiliary variables
        Vector delta_ns;
        VectorResult hs, as, ln_as;

        // The residual function
        VectorResult f;
        func(f).resize(Ns);
        grad(f).resize(Ns, Ns);

        // Initiate the Newton iterations
        for(unsigned counter = 0; counter < options$.refinement.max_iterations; ++counter)
        {
            // Assemble the equilibrium constraint vector and its partial molar derivatives
            hs = balance$(state);
            grad(hs) = cols(idx_equilibrium_stable_species, grad(hs));

            // Assemble the vector of activities of the stable species and their partial molar derivatives
            as = state.activities();
            func(as) = rows(idx_stable_species, func(as));
            grad(as) = submatrix(idx_stable_species, idx_stable_species, grad(as));

            // Assemble the vector of log of the activities of the stable species and their partial molar derivatives
            func(ln_as) = func(as).array().log();
            grad(ln_as) = func(as).asDiagonal().inverse() * grad(as);

            // Assemble the residual function
            func(f).segment(00, Nc) = func(hs);
            func(f).segment(Nc, Ms) = vs * func(ln_as) - ln_Ks;

            // Assemble the Jacobian of the residual function
            grad(f).topRows(Nc)    = grad(hs);
            grad(f).bottomRows(Ms) = vs * grad(ln_as);

            // Calculate the Newton step
            delta_ns = grad(f).lu().solve(-func(f));

            // Compute the largest boundary step that does not lead to negative values
            const double alpha = calculateLargestBoundaryStep(ns, delta_ns);

            // Compute the new molar composition of the stable equilibrium species
            ns += alpha * delta_ns;

            // Update the molar composition of the stable equilibrium species
            state.setComposition(idx_stable_species, ns);

            // Update the result variable
            ++result.num_constraint_evals;
            ++result.num_iterations;
            ++result.num_objective_evals;

            // Calculate the current residual of the calculation
            const double residual = func(f).lpNorm<Infinity>();

            // Check convergence of the refinement algorithm
            if(residual < options$.refinement.tolerance)
            {
                result.converged = true;
                break;
            }
        }

        VectorResult a, ae, au, ln_au;

        a = state.activities();

        func(ae) = partitioning$.equilibriumRows(func(a));
        grad(ae) = partitioning$.equilibriumRowsCols(grad(a));

        func(as) = rows(idx_stable_species, func(a));
        grad(as) = submatrix(idx_stable_species, idx_stable_species, grad(a));

        func(au) = rows(idx_unstable_species, func(a));
        grad(au) = submatrix(idx_unstable_species, idx_unstable_species, grad(a));

        func(ln_as) = func(as).array().log();
        grad(ln_as) = func(as).asDiagonal().inverse() * grad(as);

        func(ln_au) = func(au).array().log();
        grad(ln_au) = func(au).asDiagonal().inverse() * grad(au);

        const Vector us = us0 + RT*func(ln_as);
        const Vector uu = uu0 + RT*func(ln_au);

        Vector ye = Hs.transpose().fullPivLu().solve(us);

        Vector zu = uu - Hu.transpose() * ye;

        const unsigned Ne = partitioning$.numEquilibriumSpecies();

        Vector ze = zeros(Ne);
        setRows(idx_equilibrium_unstable_species, zu, ze);

        lagrange.setMultipliers(state, ye, ze);

        Vector ln_ae = func(ae).array().log();
        Vector ue0 = partitioning$.equilibriumRows(u0);

        return result;
    }

    auto refine(ChemicalState& state, EquilibriumLagrange& lagrange, const EquilibriumConstraints& constraints) -> EquilibriumResult
    {
        EquilibriumResult result;

        // The temperature of the chemical system
        const double T = state.temperature();

        // The pressure of the chemical system
        const double P = state.pressure();

        // The product of temperature and the universal gas constant
        const double RT = universalGasConstant * T;

        // The indices of the stable equilibrium species (indices relative to the set of all species)
        const Indices& idx_stable_species = lagrange.idxStableSpecies();

        // The indices of the unstable equilibrium species (indices relative to the set of all species)
        const Indices& idx_unstable_species = lagrange.idxUnstableSpecies();

        // The indices of the stable equilibrium species (indices relative to the set of equilibrium species)
        const Indices& idx_equilibrium_stable_species = lagrange.idxEquilibriumStableSpecies();

        // The number of equilibrium constraints
        const unsigned Nc = constraints.numConstraints();

        // The number of stable equilibrium species
        const unsigned Ns = idx_stable_species.size();

        // Set the molar composition of the unstable species to zero
        state.setComposition(idx_unstable_species, zeros(idx_unstable_species.size()));

        // Calculate the standard chemical potential of the stable equilibrium species
        const Vector us0 = rows(idx_stable_species, system$.chemicalPotentials(T, P));

        // The scaling factor used in the calculation
        const double scaling = RT;

        // The molar composition of the stable equilibrium species
        Vector ns = rows(idx_stable_species, state.composition());

        // Auxiliary variables
        Vector delta, delta_ns;
        VectorResult hs, as, ln_as, us;

        // The residual function
        VectorResult f;
        func(f).resize(Ns + Nc);
        grad(f).resize(Ns + Nc, Ns + Nc);

        // The Lagrange multipliers y and z
        Vector ys = lagrange.lagrangeY() / scaling;
        Vector zs = lagrange.lagrangeZ();

        // Initiate the Newton iterations
        for(unsigned counter = 0; counter < options$.refinement.max_iterations; ++counter)
        {
            // Assemble the equilibrium constraint vector and its partial molar derivatives
            hs = constraints(state);
            grad(hs) = cols(idx_equilibrium_stable_species, grad(hs));

            // Assemble the vector of activities of the stable species and their partial molar derivatives
            as = state.activities();
            func(as) = rows(idx_stable_species, func(as));
            grad(as) = submatrix(idx_stable_species, idx_stable_species, grad(as));

            // Assemble the vector of log of the activities of the stable species and their partial molar derivatives
            func(ln_as) = func(as).array().log();
            grad(ln_as) = func(as).asDiagonal().inverse() * grad(as);

            // Assemble the vector of chemical potentials of the stable species and their partial molar derivatives
            func(us) = (us0 + RT * func(ln_as))/scaling;
            grad(us) = (RT * grad(ln_as))/scaling;

            // Assemble the residual function
            func(f).segment(00, Ns) = func(us) + grad(hs).transpose() * ys;
            func(f).segment(Ns, Nc) = func(hs);

            // Assemble the Jacobian of the residual function
            grad(f).    topLeftCorner(Ns, Ns) = grad(us);
            grad(f).   topRightCorner(Ns, Nc) = grad(hs).transpose();
            grad(f). bottomLeftCorner(Nc, Ns) = grad(hs);
            grad(f).bottomRightCorner(Nc, Nc) = zeros(Nc, Nc);

            // Calculate the Newton step
            delta = grad(f).lu().solve(-func(f));

            // Extract the Newton step components of the stable species
            delta_ns = delta.segment(0, Ns);

            // Compute the largest boundary step that does not lead to negative values
            const double alpha = calculateLargestBoundaryStep(ns, delta_ns);

            // Compute the new molar composition of the stable equilibrium species
            ns += alpha * delta_ns;
            ys += alpha * delta.segment(Ns, Nc);

            // Update the molar composition of the stable equilibrium species
            state.setComposition(idx_stable_species, ns);

            // Update the result variable
            ++result.num_constraint_evals;
            ++result.num_iterations;
            ++result.num_objective_evals;

            // Calculate the current residual of the calculation
            const double residual = func(f).lpNorm<Infinity>();

            // Check convergence of the refinement algorithm
            if(residual < options$.refinement.tolerance)
            {
                result.converged = true;
                break;
            }
        }


        // Set the refined Lagrange multipliers y
        lagrange.setMultipliers(state, ys * scaling, zs);

        return result;
    }

    auto solve(ChemicalState& state, EquilibriumLagrange& lagrange, const EquilibriumConstraints& constraints) -> EquilibriumResult
    {
        EquilibriumResult result1 = minimise(state, lagrange, constraints);
        EquilibriumResult result2 = refine(state, lagrange);
        return result1 + result2;
    }

    auto solve(ChemicalState& state, EquilibriumLagrange& lagrange, const Vector& be) -> EquilibriumResult
    {
        EquilibriumConstraints constraints = massBalance(system$, partitioning$, be);
        return solve(state, lagrange, constraints);
    }
};

EquilibriumSolver::EquilibriumSolver(const ChemicalSystem& system)
: pimpl(new Impl(system, Partitioning(system)))
{}

EquilibriumSolver::EquilibriumSolver(const ChemicalSystem& system, const Partitioning& partitioning)
: pimpl(new Impl(system, partitioning))
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
    pimpl->setOptions(options);
}

auto EquilibriumSolver::setParams(const EquilibriumParams& params) -> void
{
    pimpl->setParams(params);
}

auto EquilibriumSolver::setScaling(const Optima::Scaling& scaling) -> void
{
    pimpl->setScaling(scaling);
}

auto EquilibriumSolver::setScaling(const ChemicalState& state) -> void
{
    pimpl->setScaling(state);
}

auto EquilibriumSolver::minimise(ChemicalState& state, EquilibriumLagrange& lagrange, const EquilibriumConstraints& constraints) -> EquilibriumResult
{
    return pimpl->minimise(state, lagrange, constraints);
}

auto EquilibriumSolver::minimise(ChemicalState& state, EquilibriumLagrange& lagrange, const Vector& be) -> EquilibriumResult
{
    return pimpl->minimise(state, lagrange, be);
}

auto EquilibriumSolver::refine(ChemicalState& state, EquilibriumLagrange& lagrange) -> EquilibriumResult
{
    return pimpl->refine(state, lagrange);
}

auto EquilibriumSolver::refine(ChemicalState& state, EquilibriumLagrange& lagrange, const EquilibriumConstraints& constraints) -> EquilibriumResult
{
    return pimpl->refine(state, lagrange, constraints);
}

auto EquilibriumSolver::solve(ChemicalState& state, EquilibriumLagrange& lagrange, const EquilibriumConstraints& constraints) -> EquilibriumResult
{
    return pimpl->solve(state, lagrange, constraints);
}

auto EquilibriumSolver::solve(ChemicalState& state, EquilibriumLagrange& lagrange, const Vector& be) -> EquilibriumResult
{
    return pimpl->solve(state, lagrange, be);
}

} // namespace Reaktor
