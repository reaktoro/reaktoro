// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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
#include <Reaktoro/Common/ChemicalVector.hpp>
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/ConvertUtils.hpp>
#include <Reaktoro/Common/Exception.hpp>
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

    /// The solver for the optimisation calculations
    OptimumSolver solver;

    /// The molar amounts of the species
    Vector n;

    /// The dual potentials of the elements
    Vector y;

    /// The dual potentials of the species
    Vector z;

    /// The molar amounts of the elements in the equilibrium partition
    Vector be;

    /// The chemical potentials of the species
    ChemicalVector u;

    /// The chemical potentials of the equilibrium species
    ChemicalVector ue;

    /// The molar fractions of the equilibrium species
    ChemicalVector xe;

    /// The optimisation problem
    OptimumProblem optimum_problem;

    /// The state of the optimisation calculation
    OptimumState optimum_state;

    /// The options for the optimisation calculation
    OptimumOptions optimum_options;

    /// The indices of the equilibrium species
    Indices ies;

    /// The indices of the equilibrium elements
    Indices iee;

    /// The number of species and elements in the system
    unsigned N, E;

    /// The number of species and elements in the equilibrium partition
    unsigned Ne, Ee;

    /// The formula matrix of the species in the system
    Matrix A;

    /// The formula matrix of the species in the equilibrium partition
    Matrix Ae;

    /// Construct a default Impl instance
    Impl()
    {}

    /// Construct a Impl instance
    Impl(const ChemicalSystem& system)
    : system(system)
    {
        // Initialize the formula matrix
        A = system.formulaMatrix();

        // Initialize the number of species and elements in the system
        N = system.numSpecies();
        E = system.numElements();

        // Set the default partition as all species are in equilibrium
        setPartition(Partition(system));
    }

    /// Set the partition of the chemical system
    auto setPartition(const Partition& partition_) -> void
    {
        // Set the partition of the chemical system
        partition = partition_;

        // Initialize the number of species and elements in the equilibrium partition
        Ne = partition.numEquilibriumSpecies();
        Ee = partition.numEquilibriumElements();

        // Initialize the formula matrix of the equilibrium species
        Ae = partition.formulaMatrixEquilibriumPartition();

        // Initialize the indices of the equilibrium species and elements
        ies = partition.indicesEquilibriumSpecies();
        iee = partition.indicesEquilibriumElements();
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

        // Set a non-zero value to the maximum denominator in the regularized balance matrix
        if(optimum_options.regularization.max_denominator == 0)
            optimum_options.regularization.max_denominator = 1e6;
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

        // The thermodynamic properties of the chemical system
        ChemicalProperties properties;

        // The result of the objective evaluation
        ObjectiveResult res;

        // The normalized standard Gibbs energies of the species at (T,P)
        ThermoVector G0 = system.properties(T, P).
            standardPartialMolarGibbsEnergies()/RT;

        // The Gibbs energy function to be minimized
        optimum_problem.objective = [=](const Vector& ne) mutable
        {
            // Set the molar amounts of the species
            rows(n, ies) = ne;

            // Calculate the thermodynamic properties of the chemical system
            properties = system.properties(T, P, n);

            // Set the scaled chemical potentials of the species
            u = G0 + properties.lnActivities();

            // Set the scaled chemical potentials of the equilibrium species
            ue = u.rows(ies, ies);

            // Set the molar fractions of the equilibrium species
            xe = properties.molarFractions().rows(ies, ies);

            // Set the objective result
            res.val = dot(ne, ue.val);
            res.grad = ue.val;

            // Set the Hessian of the objective function
            switch(options.hessian)
            {
            case GibbsHessian::Exact:
                res.hessian.mode = Hessian::Dense;
                res.hessian.dense = ue.ddn;
                break;
            case GibbsHessian::ExactDiagonal:
                res.hessian.mode = Hessian::Diagonal;
                res.hessian.diagonal = diagonal(ue.ddn);
                break;
            case GibbsHessian::Approximation:
                res.hessian.mode = Hessian::Dense;
                res.hessian.dense = diag(inv(xe.val)) * xe.ddn;
                break;
            case GibbsHessian::ApproximationDiagonal:
                res.hessian.mode = Hessian::Diagonal;
                res.hessian.diagonal = diagonal(xe.ddn)/xe.val;
                break;
            }

            return res;
        };

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
        y = state.elementDualPotentials()/RT;

        // Set the normalized dual potentials of the species
        z = state.speciesDualPotentials()/RT;

        // Initialize the optimum state
        rows(n, ies).to(optimum_state.x);
        rows(y, iee).to(optimum_state.y);
        rows(z, ies).to(optimum_state.z);
    }

    /// Initialize the chemical state from a optimum state
    auto updateChemicalState(ChemicalState& state) -> void
    {
        // The temperature and the RT factor
        const double T  = state.temperature();
        const double RT = universalGasConstant*T;

        // Update the molar amounts of the equilibrium species
        rows(n, ies) = optimum_state.x;

        // Update the dual potentials of the species and elements (in units of J/mol)
        z = zeros(N); rows(z, ies) = optimum_state.z * RT;
        y = zeros(E); rows(y, iee) = optimum_state.y * RT;

        // Update the chemical state
        state.setSpeciesAmounts(n);
        state.setElementDualPotentials(y);
        state.setSpeciesDualPotentials(z);
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
        y = state.elementDualPotentials();
        z = state.speciesDualPotentials();

        // Calculate the standard thermodynamic properties of the system
        ChemicalProperties props = system.properties(T, P, n);

        // Get the standard Gibbs energies of the equilibrium species
        const Vector ge0 = rows(props.standardPartialMolarGibbsEnergies().val, ies);

        // Get the ln activity constants of the equilibrium species
        const Vector ln_ce = rows(props.lnActivityConstants().val, ies);

        // Define the optimisation problem
        OptimumProblem optimum_problem;
        optimum_problem.c = ge0/RT + ln_ce;
        optimum_problem.A = Ae;
        optimum_problem.b = be;
        optimum_problem.l = zeros(Ne);
        optimum_problem.u = ones(Ne) * inf;

        // Initialize the optimum state
        OptimumState optimum_state;

        // The result of the linear programming calculation
        EquilibriumResult result;

        // Initialize the solver with a simplex method
        OptimumSolver solver(OptimumMethod::Simplex);

        // Do not allow echelonization for simplex calculations.
        optimum_options.regularization.echelonize = false;

        // Solve the linear programming problem
        result.optimum = solver.solve(optimum_problem, optimum_state, optimum_options);

        // Update the molar amounts of the equilibrium species
        rows(n, ies) = optimum_state.x;

        // Update the dual potentials of the species and elements (in units of J/mol)
        z = zeros(N); rows(z, ies) = optimum_state.z * RT;
        y = zeros(E); rows(y, iee) = optimum_state.y * RT;

        // Update the chemical state
        state.setSpeciesAmounts(n);
        state.setElementDualPotentials(y);
        state.setSpeciesDualPotentials(z);

        return result;
    }

    /// Find an initial guess for an equilibrium problem.
    auto initialguess(ChemicalState& state, double T, double P, Vector be) -> EquilibriumResult
    {
        // Solve the linear programming problem to obtain an approximation
        auto result = approximate(state, T, P, be);

        // Auxiliary variables
        const double RT = universalGasConstant*T;

        // Replace zero amounts by a positive small amount
        for(Index i : ies)
        {
            // This is an extremely part of settint initial guess.
            // The value 1e-12 is to clean any entry that is
            // contaminated with round-off errors. The value 1e-10
            // is used as an initial amount for species with zero amounts.
            // It is neither an extremely very small value that
            // causes some primal variables to get prematurely
            // trapped on the bounds, nor a too big value that
            // spoils the mass balance residuals obtained from
            // the simplex calculation. The RT factor is to put
            // the Lagrange multipliers z with J/mol scale.
            // Be very carefull in changing these values!
            n[i] = (n[i] > 1e-12) ? n[i] : 1e-10;
            z[i] = (z[i] > 1e-12) ? z[i] : 1e-10 * RT;
        }

        // Update the chemical state
        state.setSpeciesAmounts(n);
        state.setElementDualPotentials(y);
        state.setSpeciesDualPotentials(z);

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

    /// Solve the equilibrium problem
    auto solve(ChemicalState& state, double T, double P, const Vector& be) -> EquilibriumResult
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
    auto solve(ChemicalState& state, double T, double P, const double* b) -> EquilibriumResult
    {
        // Set the molar amounts of the elements
        be = Vector::Map(b, Ee);

        // Set temperature and pressure of the chemical state
        state.setTemperature(T);
        state.setPressure(P);

        // Check if a simplex cold-start approximation must be performed
        if(coldstart(state))
            initialguess(state, T, P, be);

        // The result of the equilibrium calculation
        EquilibriumResult result;

        // Update the optimum options
        updateOptimumOptions();

        // Update the optimum problem
        updateOptimumProblem(state);

        // Update the optimum state
        updateOptimumState(state);

        // Set the method for the optimisation calculation
        solver.setMethod(options.method);

        // Solve the optimisation problem
        result.optimum += solver.solve(optimum_problem, optimum_state, optimum_options);

        // Update the chemical state from the optimum state
        updateChemicalState(state);

        return result;
    }

    /// Return the sensitivity of the equilibrium state.
    auto sensitivity() -> EquilibriumSensitivity
    {
        Vector zerosEe = zeros(Ee);
        Vector zerosNe = zeros(Ne);
        Vector unitjEe = zeros(Ee);

        EquilibriumSensitivity sensitivity;
        sensitivity.dnedT.resize(Ne);
        sensitivity.dnedP.resize(Ne);
        sensitivity.dnedbe.resize(Ne, Ee);

        sensitivity.dnedT = solver.dxdp(ue.ddt, zerosEe);
        sensitivity.dnedP = solver.dxdp(ue.ddp, zerosEe);
        for(Index j = 0; j < Ee; ++j)
        {
            unitjEe = unit(Ee, j);
            sensitivity.dnedbe.col(j) = solver.dxdp(zerosNe, unitjEe);
        }

        return sensitivity;
    }
};

EquilibriumSolver::EquilibriumSolver()
: pimpl(new Impl())
{}

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

auto EquilibriumSolver::approximate(ChemicalState& state, double T, double P, const Vector& be) -> EquilibriumResult
{
    return pimpl->approximate(state, T, P, be);
}

auto EquilibriumSolver::solve(ChemicalState& state, double T, double P, const Vector& be) -> EquilibriumResult
{
    return pimpl->solve(state, T, P, be);
}

auto EquilibriumSolver::solve(ChemicalState& state, double T, double P, const double* be) -> EquilibriumResult
{
    return pimpl->solve(state, T, P, be);
}

auto EquilibriumSolver::sensitivity() -> EquilibriumSensitivity
{
    return pimpl->sensitivity();
}

} // namespace Reaktoro
