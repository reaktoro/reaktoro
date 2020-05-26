// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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

// Optima includes
#include <Optima/Options.hpp>
#include <Optima/Problem.hpp>
#include <Optima/Result.hpp>
#include <Optima/Solver.hpp>
#include <Optima/State.hpp>

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConstraints.hpp>
#include <Reaktoro/Equilibrium/EquilibriumDims.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>

namespace Reaktoro {
namespace detail {


// auto createOptimaProblem(const EquilibriumConstraints& constraints) -> Optima::Problem
// {

// }

// auto createOptimaProblem(const EquilibriumConstraints& constraints) -> Optima::Problem
// {

// }

// auto updateOptimaProblem(const EquilibriumConstraints& constraints, Optima::Problem& problem)
// {

// }




/// Return the Optima::Dims object with dimensions of the optimization problem.
auto dims(const EquilibriumProblem& problem) -> Optima::Dims
{
    Optima::Dims dims;
    dims.x = problem.dims().Nx;
    dims.be = problem.dims().Nc;
    return dims;
}

} // namespace detail

struct EquilibriumSolver::Impl
{
    /// The chemical system associated with this equilibrium solver
    ChemicalSystem system;

    /// The equilibrium constraints associated with this equilibrium solver
    EquilibriumConstraints constraints;

    /// The equilibrium problem with conservation matrix and objective function
    EquilibriumProblem eqproblem;

    /// The options of the equilibrium solver
    EquilibriumOptions eqoptions;

    /// The optimization problem corresponding to the constrained equilibrium calculation.
    Optima::Problem optproblem;

    /// The optimization state of the calculation
    Optima::State optstate;

    /// The solver for the optimization calculations
    Optima::Solver optsolver;

    /// The result of the optimization calculations
    Optima::Result optresult;

    /// The options of the optimization solver
    Optima::Options optoptions;

    Index Ne  = 0; ///< The number of chemical elements and electric charge in the chemical system.
    Index Nn  = 0; ///< The number of species in the chemical system.
    Index Npe = 0; ///< The number of equation constraints among the functional constraints.
    Index Npp = 0; ///< The number of property preservation constraints among the functional constraints.
    Index Np  = 0; ///< The number of functional constraints (Np = Npe + Npp).
    Index Nq  = 0; ///< The number of chemical potential constraints.
    Index Nir = 0; ///< The number of reactions prevented from reacting during the equilibrium calculation.
    Index Nc  = 0; ///< The number of introduced control variables (must be equal to Np)
    Index Nx  = 0; ///< The number of variables (Nx = Nn + Np + Nq).
    Index Nb  = 0; ///< The number of components (Nb = Ne + Nir).

    /// Construct a Impl instance with given EquilibriumConstraints object
    Impl(const EquilibriumConstraints& constraints)
    : system(constraints.system()),
      constraints(constraints),
      eqproblem(constraints),
      optproblem(detail::dims(eqproblem)),
      optstate(detail::dims(eqproblem)),
      optsolver(optproblem)
    {
        // Initialize the auxiliary number variables
        Ne  = system.elements().size() + 1;
        Nn  = system.species().size();
        Npe = constraints.data().econstraints.size();
        Npp = constraints.data().pconstraints.size();
        Np  = Npe + Npp;
        Nq  = constraints.data().uconstraints.size();
        Nir = constraints.data().restrictions.reactions_cannot_react.size();
        Nc  = constraints.data().controls.size();
        Nx  = Nn + Np + Nq;
        Nb  = Ne + Nir;
    }

    /// Update the equilibrium constraints before a new calculation
    auto updateConstraints(const EquilibriumConstraints& constraints)
    {
        // Update the equilibrium problem with updated equilibrium constraints
        eqproblem.update(constraints);
    }

    /// Update the optimization problem before a new equilibrium calculation.
    auto updateOptProblem(const ChemicalProps& props0)
    {
        // Create the objective function corresponding to the constrained equilibrium calculation
        EquilibriumObjective obj = eqproblem.objective(props0);

        // Update the objective function in the Optima::Problem object
        optproblem.f = [=](VectorXdConstRef x, Optima::ObjectiveResult& res)
        {
            if(res.requires.f) res.f = obj.f(x);
            if(res.requires.g) obj.g(x, res.g);
            if(res.requires.H) obj.H(x, res.H);
        };
    }

    /// Update the initial state variables before the new equilibrium calculation.
    auto updateOptState(const ChemicalProps& props0)
    {
        auto n = optstate.x.head(Nn);
        auto p = optstate.x.segment(Nn, Np);
        auto q = optstate.x.tail(Nq);

        // n = props0.speciesAmounts();
        // p =
    }

    /// Return the objective function in Optima's format.
    auto objectiveFunction(const ChemicalProps& props0) const -> Optima::ObjectiveFunction
    {
        EquilibriumObjective obj = eqproblem.objective(props0);

        return [=](VectorXdConstRef x, Optima::ObjectiveResult& res)
        {
            if(res.requires.f) res.f = obj.f(x);
            if(res.requires.g) obj.g(x, res.g);
            if(res.requires.H) obj.H(x, res.H);
        };
    }

    /// Solve an equilibrium problem with given chemical state in disequilibrium.
    auto solve(ChemicalState& state) -> EquilibriumResult
    {
        EquilibriumResult eqresult;

        ChemicalProps props0 = state.props();

        optproblem.f = objectiveFunction(props0);

        // updateOptOptions(options);
        // updateOptProblem(props0);
        // updateOptState(props0);

        optresult = optsolver.solve(optstate, optproblem);

        return eqresult;
    }

//     /// Update the OptimumOptions instance with given EquilibriumOptions instance
//     auto updateOptimumOptions() -> void
//     {
//         // Initialize the options for the optimization calculation
//         optimum_options = options.optimum;

//         // Set the parameters of the optimization algorithms that control how small can be the amount of a species
//         optimum_options.ipaction.mu = options.epsilon;
//         optimum_options.ipnewton.mu = options.epsilon;
//         optimum_options.ipopt.mu.push_back(options.epsilon);
//         optimum_options.ipactive.epsilon = options.epsilon;

//         // Initialize the names of the primal and dual variables
//         if(options.optimum.output.active)
//         {
//             // Use `n` instead of `x` to name the variables
//             optimum_options.output.xprefix = "n";

//             // Define some auxiliary references to the variables names
//             auto& xnames = optimum_options.output.xnames;
//             auto& ynames = optimum_options.output.ynames;
//             auto& znames = optimum_options.output.znames;

//             // Initialize the names of the primal variables `n`
//             for(Index i : ies)
//                 xnames.push_back(system.species(i).name());

//             // Initialize the names of the dual variables `y`
//             for(Index i : iee)
//                 ynames.push_back(system.element(i).name());

//             // Initialize the names of the dual variables `z`
//             znames = xnames;
//         }
//     }

//     /// Update the OptimumProblem instance with given EquilibriumProblem and ChemicalState instances
//     auto updateOptimumProblem(const ChemicalState& state) -> void
//     {
//         // The temperature and pressure of the equilibrium calculation
//         const auto T  = state.temperature();
//         const auto P  = state.pressure();
//         const auto RT = universalGasConstant*T;

//         // Set the molar amounts of the species
//         n = state.speciesAmounts();

//         // Update the standard thermodynamic properties of the chemical system
//         props.update(T, P, n);

//         // Update the normalized standard Gibbs energies of the species
//         u0 = props.standardGibbsEnergies()/RT;

//         // The result of the objective evaluation
//         ObjectiveResult res;

//         // The Gibbs energy function to be minimized
//         optimum_problem.objective = [=](ArrayXrConstRef ne) mutable
//         {
//             // Set the molar amounts of the species
//             n(ies) = ne;

//             // Update the chemical properties of the chemical system
//             props.update(T, P, n);

//             // Set the scaled chemical potentials of the species
//             u = u0 + props.lnActivities();

//             // Set the scaled chemical potentials of the equilibrium species
//             ue = u(ies);

//             // Set the mole fractions of the equilibrium species
//             xe = props.moleFractions()(ies);

//             // Set the objective result
//             res.val = (ne * ue).sum();
//             res.grad = ue;

//             // Set the Hessian of the objective function
//             switch(options.hessian)
//             {
//             case GibbsHessian::Exact:
//                 // res.hessian.mode = Hessian::Dense;
//                 // res.hessian.dense = ue.ddn;
//                 // break;
//             case GibbsHessian::ExactDiagonal:
//                 // res.hessian.mode = Hessian::Diagonal;
//                 // res.hessian.diagonal = diagonal(ue.ddn);
//                 // break;
//             case GibbsHessian::Approximation:
//                 // res.hessian.mode = Hessian::Dense;
//                 // res.hessian.dense = diag(inv(xe.val)) * xe.ddn;
//                 // break;
//             case GibbsHessian::ApproximationDiagonal:
//                 res.hessian.mode = Hessian::Diagonal;
//                 res.hessian.diagonal = 1.0 / ne.array();
//                 break;
//             }

//             return res;
//         };

//         optimum_problem.c.resize(0);
//         optimum_problem.n = Ne;
//         optimum_problem.A = Ae;
//         optimum_problem.b = be;
//         optimum_problem.l.setConstant(Ne, options.epsilon);
//     }

//     /// Initialize the optimum state from a chemical state
//     auto updateOptimumState(const ChemicalState& state) -> void
//     {
//         // The temperature and the RT factor
//         const double T  = state.temperature();
//         const double RT = universalGasConstant*T;

//         // Set the molar amounts of the species
//         n = state.speciesAmounts();

//         // Set the normalized dual potentials of the elements
//         y = state.elementDualPotentials()/RT;

//         // Set the normalized dual potentials of the species
//         z = state.speciesDualPotentials()/RT;

//         // Initialize the optimum state
//         optimum_state.x = n(ies);
//         optimum_state.y = y(iee);
//         optimum_state.z = z(ies);
//     }

//     /// Initialize the chemical state from a optimum state
//     auto updateChemicalState(ChemicalState& state) -> void
//     {
//         // The temperature and the RT factor
//         const double T  = state.temperature();
//         const double RT = universalGasConstant*T;

//         // Update the molar amounts of the equilibrium species
//         n(ies) = optimum_state.x;

//         // Update the normalized chemical potentials of the inert species
//         ui = u(iis);

//         // Update the normalized dual potentials of the elements
//         y = zeros(E); y(iee) = optimum_state.y;

//         // Update the normalized dual potentials of the equilibrium and inert species
//         z(ies) = optimum_state.z;
//         z(iis) = ui - tr(Ai) * y;

//         // Scale the normalized dual potentials of elements and species to units of J/mol
//         y *= RT;
//         z *= RT;

//         // Update the chemical state
//         state.setSpeciesAmounts(n);
//         state.setElementDualPotentials(y);
//         state.setSpeciesDualPotentials(z);
//     }

//     /// Find a feasible approximation for an equilibrium problem with all elements present on chemical system.
//     auto approximate_with_all_element_amounts(ChemicalState& state, double T, double P, ArrayXrConstRef b) -> EquilibriumResult
//     {
//         ArrayXr be = b(iee);
//         return approximate(state, T, P, be);
//     }

//     /// Find a feasible approximation for an equilibrium problem.
//     auto approximate(ChemicalState& state, double T, double P, ArrayXr be) -> EquilibriumResult
//     {
//         // Check the dimension of the vector `be`
//         Assert(unsigned(be.rows()) == Ee,
//             "Cannot proceed with method EquilibriumSolver::approximate.",
//             "The dimension of the given vector of molar amounts of the "
//             "elements does not match the number of elements in the "
//             "equilibrium partition.");

//         // Set temperature and pressure of the chemical state
//         state.setTemperature(T);
//         state.setPressure(P);

//         // Auxiliary variables
//         const double RT = universalGasConstant*T;
//         const double inf = std::numeric_limits<double>::infinity();

//         // Update the internal state of n, y, z
//         n = state.speciesAmounts();
//         y = state.elementDualPotentials();
//         z = state.speciesDualPotentials();

//         // Update the standard thermodynamic properties of the system
//         props.update(T, P);

//         // Get the standard Gibbs energies of the equilibrium species
//         const ArrayXr ge0 = props.standardGibbsEnergies()(ies);

//         // Get the ln activity constants of the equilibrium species
//         const ArrayXr ln_ce = props.lnActivityConstants()(ies);

//         // Define the optimization problem
//         OptimumProblem optimum_problem;
//         optimum_problem.n = Ne;
//         optimum_problem.c = ge0/RT + ln_ce;
//         optimum_problem.A = Ae;
//         optimum_problem.b = be;
//         optimum_problem.l = zeros(Ne);
//         optimum_problem.u = ones(Ne) * inf;

//         // Initialize the optimum state
//         OptimumState optimum_state;

//         // The result of the linear programming calculation
//         EquilibriumResult result;

//         // Update the optimum options
//         updateOptimumOptions();

//         // Set the method for the optimization calculation
//         solver.setMethod(options.method);

//         // Solve the linear programming problem
//         result.optimum = solver.solve(optimum_problem, optimum_state, optimum_options);

//         // Update the molar amounts of the equilibrium species
//         n(ies) = optimum_state.x;

//         // Update the dual potentials of the species and elements (in units of J/mol)
//         z = zeros(N); z(ies) = optimum_state.z * RT;
//         y = zeros(E); y(iee) = optimum_state.y * RT;

//         // Update the chemical state
//         state.setSpeciesAmounts(n);
//         state.setElementDualPotentials(y);
//         state.setSpeciesDualPotentials(z);

//         return result;
//     }

//     /// Find an initial guess for an equilibrium problem.
//     auto initialguess(ChemicalState& state, double T, double P, ArrayXr be) -> EquilibriumResult
//     {
//         // Solve the linear programming problem to obtain an approximation
//         auto result = approximate(state, T, P, be);

//         // Check the approximate calculation was successful
//         Assert(result.optimum.succeeded,
//             "Cannot proceed with the equilibrium calculation.",
//             "The calculation of initial guess failed, most "
//             "probably because no feasible solution exists.");

//         // Preserve the values of n that are greater than z. For all others, set n to sqrt(epsilon)
//         n = (n.array() > z.array()).select(n, std::sqrt(options.epsilon));

//         // Set z to zero, which will later be set to epsilon / n.
//         z.fill(0.0);

//         // Update the chemical state
//         state.setSpeciesAmounts(n);
//         state.setSpeciesDualPotentials(z);

//         return result;
//     }

//     /// Return true if cold-start is needed.
//     auto coldstart(const ChemicalState& state) -> bool
//     {
//         // Check if all equilibrium species have zero amounts
//         bool zero = true;
//         for(Index i : ies)
//             if(state.speciesAmount(i) > 0)
//                 { zero = false; break; }
//         return zero || !options.warmstart;
//     }

//     /// Solve the equilibrium problem, passing all elements that has on chemical system
//     auto solve_with_all_element_amounts(ChemicalState& state, double T, double P, ArrayXrConstRef b) -> EquilibriumResult
//     {
//         ArrayXr be = b(iee);
//         return solve(state, T, P, be);
//     }

//     /// Solve the equilibrium problem
//     auto solve(ChemicalState& state, double T, double P, ArrayXrConstRef be) -> EquilibriumResult
//     {
//         // Check the dimension of the vector `be`
//         Assert(be.size() == static_cast<int>(Ee),
//             "Cannot proceed with method EquilibriumSolver::solve.",
//             "The dimension of the given vector of molar amounts of the "
//             "elements does not match the number of elements in the "
//             "equilibrium partition.");
//         return solve(state, T, P, be.data());
//     }

//     /// Solve the equilibrium problem
//     auto solve(ChemicalState& state, double T, double P, const double* _be) -> EquilibriumResult
//     {
//         // Set the molar amounts of the elements
//         be = ArrayXr::Map(_be, Ee);

//         // Set temperature and pressure of the chemical state
//         state.setTemperature(T);
//         state.setPressure(P);

//         // Check if a simplex cold-start approximation must be performed
//         if(coldstart(state))
//             initialguess(state, T, P, be);

//         // The result of the equilibrium calculation
//         EquilibriumResult result;

//         // Update the optimum options
//         updateOptimumOptions();

//         // Update the optimum problem
//         updateOptimumProblem(state);

//         // Update the optimum state
//         updateOptimumState(state);

//         // Set the method for the optimization calculation
//         solver.setMethod(options.method);

//         // Set the maximum number of iterations in each optimization pass
//         optimum_options.max_iterations = 10;

//         // Start the several opmization passes (stop if convergence attained)
//         auto counter = 0;
//         while(counter < options.optimum.max_iterations)
//         {
//             // Solve the optimization problem
//             result.optimum += solver.solve(optimum_problem, optimum_state, optimum_options);

//             // Exit this loop if last solve succeeded
//             if(result.optimum.succeeded)
//                 break;

//             counter += optimum_options.max_iterations;
//         }

//         // Update the chemical state from the optimum state
//         updateChemicalState(state);

//         return result;
//     }

//     /// Return the sensitivity of the equilibrium state.
//     auto sensitivity() -> const EquilibriumSensitivity&
//     {
//         zerosEe = zeros(Ee);
//         zerosNe = zeros(Ne);
//         unitjEe = zeros(Ee);

//         sensitivities.dndT = zeros(Ne);
//         sensitivities.dndP = zeros(Ne);
//         sensitivities.dndb = zeros(Ne, Ee);

//         // sensitivities.dndT = solver.dxdp(ue.ddT, zerosEe);
//         // sensitivities.dndP = solver.dxdp(ue.ddP, zerosEe);
//         for(Index j = 0; j < Ee; ++j)
//         {
//             unitjEe = unit(Ee, j);
//             sensitivities.dndb.col(j) = solver.dxdp(zerosNe, unitjEe);
//         }

//         return sensitivities;
//     }

//     /// Compute the sensitivity of the species amounts with respect to temperature.
//     auto dndT() -> ArrayXrConstRef
//     {
//         const auto& ieq_species = partition.indicesEquilibriumSpecies();
//         zerosEe = zeros(Ee);
//         sensitivities.dndT = zeros(N);
//         // sensitivities.dndT(ieq_species) = solver.dxdp(ue.ddT, zerosEe);
//         return sensitivities.dndT;
//     }

//     /// Compute the sensitivity of the species amounts with respect to pressure.
//     auto dndP() -> ArrayXrConstRef
//     {
//         const auto& ieq_species = partition.indicesEquilibriumSpecies();
//         zerosEe = zeros(Ee);
//         sensitivities.dndP = zeros(N);
//         // sensitivities.dndP(ieq_species) = solver.dxdp(ue.ddP, zerosEe);
//         return sensitivities.dndP;
//     }

//     /// Compute the sensitivity of the species amounts with respect to element amounts.
//     auto dndb() -> ArrayXrConstRef
//     {
//         const auto& ieq_species = partition.indicesEquilibriumSpecies();
//         const auto& ieq_elements = partition.indicesEquilibriumElements();
//         zerosEe = zeros(Ee);
//         zerosNe = zeros(Ne);
//         unitjEe = zeros(Ee);
//         sensitivities.dndb = zeros(Ne, Ee);
//         for(Index j : ieq_elements)
//         {
//             unitjEe = unit(Ee, j);
//             sensitivities.dndb.col(j)(ieq_species) = solver.dxdp(zerosNe, unitjEe);
//         }
//         return sensitivities.dndb;
//     }
};

EquilibriumSolver::EquilibriumSolver(const ChemicalSystem& system)
: pimpl(new Impl(EquilibriumConstraints(system)))
{}

EquilibriumSolver::EquilibriumSolver(const EquilibriumConstraints& constraints)
: pimpl(new Impl(constraints))
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
    pimpl->eqoptions = options;
}

auto EquilibriumSolver::solve(ChemicalState& state) -> EquilibriumResult
{
    return {};
}

// auto EquilibriumSolver::approximate(ChemicalState& state, double T, double P, ArrayXrConstRef be) -> EquilibriumResult
// {
//     return pimpl->approximate(state, T, P, be);
// }

// auto EquilibriumSolver::approximate(ChemicalState& state, const EquilibriumProblem& problem) -> EquilibriumResult
// {
//     setPartition(problem.partition());
//     return pimpl->approximate_with_all_element_amounts(state, problem.temperature(), problem.pressure(), problem.elementAmounts());
// }

// auto EquilibriumSolver::approximate(ChemicalState& state) -> EquilibriumResult
// {
//     return pimpl->approximate_with_all_element_amounts(state, state.temperature(), state.pressure(), state.elementAmounts());
// }

// auto EquilibriumSolver::solve(ChemicalState& state, double T, double P, ArrayXrConstRef be) -> EquilibriumResult
// {
//     return pimpl->solve(state, T, P, be);
// }

// auto EquilibriumSolver::solve(ChemicalState& state, double T, double P, const double* be) -> EquilibriumResult
// {
//     return pimpl->solve(state, T, P, be);
// }

// auto EquilibriumSolver::solve(ChemicalState& state) -> EquilibriumResult
// {
//     return pimpl->solve_with_all_element_amounts(state, state.temperature(), state.pressure(), state.elementAmounts());
// }

// auto EquilibriumSolver::solve(ChemicalState& state, const EquilibriumProblem& problem) -> EquilibriumResult
// {
//     setPartition(problem.partition());
//     return pimpl->solve_with_all_element_amounts(state, problem.temperature(), problem.pressure(), problem.elementAmounts());
// }

// auto EquilibriumSolver::props() const -> const ChemicalProps&
// {
//     return pimpl->properties;
// }

} // namespace Reaktoro
