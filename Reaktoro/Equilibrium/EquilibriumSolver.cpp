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
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/Connectivity.hpp>
#include <Reaktoro/Core/ThermoProperties.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Math/MathUtils.hpp>
#include <Reaktoro/Optimization/OptimumOptions.hpp>
#include <Reaktoro/Optimization/OptimumProblem.hpp>
#include <Reaktoro/Optimization/OptimumResult.hpp>
#include <Reaktoro/Optimization/OptimumSolver.hpp>
#include <Reaktoro/Optimization/OptimumSolverSimplex.hpp>
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

    /// The chemical potentials of the species
    ChemicalVector u;

    /// The chemical potentials of the equilibrium species
    ChemicalVector ue;

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

    /// The formula matrix of the species
    Matrix A;

    /// The formula matrix of the equilibrium species
    Matrix Ae;

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
            if(is_positive_element[index] && be[i] <= 0)
                be[i] = options.epsilon;
        }
    }

    /// Update the OptimumOptions instance with given EquilibriumOptions instance
    auto updateOptimumOptions() -> void
    {
        // Initialize the options for the optimisation calculation
        optimum_options = options.optimum;

        // Set the parameters of the optimisation algorithms that control how small can be the amount of a species
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
            for(Index i : iequilibrium_species)
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
        // The temperature and pressure of the equilibrium calculation
        const double T  = state.temperature();
        const double P  = state.pressure();
        const double RT = universalGasConstant*T;

        // Set the molar amounts of the species
        n = state.speciesAmounts();

        // The thermodynamic properties of the chemical system
        ChemicalProperties properties;

        // The result of the objective evaluation
        ObjectiveResult res;

        optimum_problem.objective = [=](const Vector& ne) mutable
        {
            // Set the molar amounts of the species
            rows(n, iequilibrium_species) = ne;

            // Calculate the thermodynamic properties of the chemical system
            properties = system.properties(T, P, n);

            // Set the scaled chemical potentials of the species
            u = properties.chemicalPotentials()/RT;

            // Set the scaled chemical potentials of the equilibrium species
            ue = u.rows(iequilibrium_species, iequilibrium_species);

            // Set the objective result
            res.val = dot(ne, ue.val);
            res.grad = ue.val;

            switch(options.hessian)
            {
            case EquilibriumHessian::Diagonal:
                res.hessian.mode = Hessian::Diagonal;
                res.hessian.diagonal = inv(ne);
                break;
            case EquilibriumHessian::SparseDiagonal:
                res.hessian.mode = Hessian::Diagonal;
                res.hessian.diagonal = inv(ne);
                for(Index i = 0; i < iequilibrium_species.size(); ++i)
                    if(system.numSpeciesInPhase(
                        system.indexPhaseWithSpecies(iequilibrium_species[i])) == 1)
                            res.hessian.diagonal[i] = 0.0;
                break;
            default:
                res.hessian.mode = Hessian::Dense;
                res.hessian.dense = ue.ddn;
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
        y = state.elementPotentials()/RT;

        // Set the normalized dual potentials of the species
        z = state.speciesPotentials()/RT;

        // Initialize the optimum state
        rows(n, iequilibrium_species).to(optimum_state.x);
        rows(y, iequilibrium_elements).to(optimum_state.y);
        rows(z, iequilibrium_species).to(optimum_state.z);
    }

    /// Initialize the chemical state from a optimum state
    auto updateChemicalState(ChemicalState& state) -> void
    {
        // The temperature and the RT factor
        const double T  = state.temperature();
        const double RT = universalGasConstant*T;

        // Update the dual potentials of the species and elements (in units of J/mol)
        z.fill(0.0); rows(z, iequilibrium_species) = optimum_state.z * RT;
        y.fill(0.0); rows(y, iequilibrium_elements) = optimum_state.y * RT;

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
        const double T = state.temperature();
        const double P = state.pressure();
        const double Pbar = P * 1e-5;
        const double RT = universalGasConstant*T;
        const double lnPbar = std::log(Pbar);
        const double inf = std::numeric_limits<double>::infinity();

        // Calculate the standard thermodynamic properties of the system
        ThermoProperties tp = system.properties(T, P);

        // Calculate the normalized standard Gibbs energies of the species
        const Vector g0 = tp.standardPartialMolarGibbsEnergies().val/RT;

        // Collect the normalized standard Gibbs energies of the equilibrium species
        const Vector ge0 = rows(g0, iequilibrium_species);

        // Initialize the connectivity of the chemical system
        Connectivity connectivity(system);

        // Compute the vector `log(ce)`, where `ce` is the vector of activity factors
        // of the equilibrium species (i.e., log(P) for species with ideal gas state
        // as standard reference)
        Vector ln_ce(Ne);
        for(Index i : iequilibrium_species)
        {
            // Get the index of the phase with current equilibrium species
            const Index iphase = connectivity.indexPhaseWithSpecies(i);

            // Check if the standard reference state of this species is ideal gas
            if(system.phase(iphase).referenceState() == PhaseReferenceState::IdealGas)
                ln_ce[i] = lnPbar;
        }

        // Define the optimisation problem
        OptimumProblem optimum_problem;
        optimum_problem.c = ge0 + ln_ce;
        optimum_problem.A = Ae;
        optimum_problem.b = be;
        optimum_problem.l = zeros(Ne);
        optimum_problem.u = ones(Ne) * inf;

        // Initialize the optimum state
        OptimumState optimum_state;

        EquilibriumResult result;

        OptimumSolverSimplex simplex;
        result.optimum = simplex.solve(optimum_problem, optimum_state);

        n = state.speciesAmounts();
        y = state.elementPotentials();
        z = state.speciesPotentials();

        // Update the chemical state from the optimum state
        rows(n, iequilibrium_species) = optimum_state.x;

        // Update the dual potentials of the species and elements
        z.fill(0.0); rows(z, iequilibrium_species) = optimum_state.z;
        y.fill(0.0); rows(y, iequilibrium_elements) = optimum_state.y;

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

        // The result of the equilibrium calculation
        EquilibriumResult result;

        // Enforce positive molar amounts for positive elements
        regularizeElementAmounts(be);

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

        return result;
    }

    /// Return the partial derivatives dn/db
    auto dndb(const ChemicalState& state) -> Matrix
    {
        const Vector ne = rows(n, iequilibrium_species);
        const Vector ze = rows(z, iequilibrium_species);

        Hessian He;
        He.mode = Hessian::Dense;
        He.dense = submatrix(u.ddn, iequilibrium_species, iequilibrium_species);

        Matrix Be;

        const Indices ili_elements = linearlyIndependentRows(Ae, Be);

        const unsigned num_li_elements = ili_elements.size();

        KktSolution sol;

        KktMatrix lhs{He, Be, ne, ze};

        KktSolver kkt;
        kkt.decompose(lhs);

        KktVector rhs;
        rhs.rx = zeros(Ne);
        rhs.rz = zeros(Ne);

        Matrix d_ne_d_be = zeros(Ne, num_li_elements);

        for(Index i = 0; i < ili_elements.size(); ++i)
        {
            rhs.ry = Vector::Unit(num_li_elements, i);
            kkt.solve(rhs, sol);
            d_ne_d_be.col(i) = sol.dx;
        }

        Matrix dndb = zeros(Ne, Ee);
        submatrix(dndb, iequilibrium_species, ili_elements) = d_ne_d_be;

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
    return pimpl->solve(state, be);
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
