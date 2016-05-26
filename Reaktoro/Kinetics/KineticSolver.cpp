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

#include "KineticSolver.hpp"

// C++ includes
#include <functional>
using namespace std::placeholders;

// Reaktoro includes
#include <Reaktoro/Common/ChemicalVector.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktoro/Kinetics/KineticOptions.hpp>
#include <Reaktoro/Kinetics/KineticProblem.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>

namespace Reaktoro {

struct KineticSolver::Impl
{
    /// The kinetically-controlled chemical reactions
    ReactionSystem reactions;

    /// The chemical system instance
    ChemicalSystem system;

    /// The partition of the species in the chemical system
    Partition partition;

    /// The options of the kinetic solver
    KineticOptions options;

    /// The equilibrium solver instance
    EquilibriumSolver equilibrium;

    /// The sensitivity of the equilibrium state
    EquilibriumSensitivity sensitivity;

    /// The ODE solver instance
    ODESolver ode;

    /// The indices of the equilibrium and kinetic species
    Indices ies, iks;

    /// The indices of the elements in the equilibrium and kinetic partition
    Indices iee, ike;

    /// The number of equilibrium and kinetic species
    unsigned Ne, Nk;

    /// The number of elements in the equilibrium and kinetic partition
    unsigned Ee, Ek;

    /// The formula matrix of the equilibrium species
    Matrix We;

    /// The stoichiometric matrix w.r.t. the equilibrium species
    Matrix Se;

    /// The stoichiometric matrix w.r.t. the kinetic species
    Matrix Sk;

    /// The coefficient matrix `A` of the chemical kinetics problem
    Matrix A;

    /// The temperature of the chemical system (in units of K)
    double T;

    /// The pressure of the chemical system (in units of Pa)
    double P;

    /// The molar composition of the equilibrium species
    Vector ne;

    /// The molar composition of the kinetic species
    Vector nk;

    /// The molar abundance of the elements in the equilibrium species
    Vector be;

    /// The combined vector of elemental molar abundance and composition of kinetic species [be nk]
    Vector benk;

    /// The chemical properties of the system
    ChemicalProperties properties;

    /// The vector with the values of the kinetic rates
    ChemicalVector rk;

    /// The partial derivatives of the reaction rates `rk` w.r.t. to `be`
    Matrix drkdbe;

    /// The partial derivatives of the reaction rates `rk` w.r.t. to `ne`
    Matrix drkdne;

    /// The partial derivatives of the reaction rates `rk` w.r.t. to `nk`
    Matrix drkdnk;

    /// The partial derivatives of the reaction rates `rk` w.r.t. to `u = [be nk]`
    Matrix drkdu;

    Impl()
    {}

    Impl(const ReactionSystem& reactions)
    : reactions(reactions), system(reactions.system()), equilibrium(system)
    {
        setPartition(Partition(system));
    }

    auto setOptions(const KineticOptions& options_) -> void
    {
        // Initialise the options of the kinetic solver
        options = options_;
    }

    auto setPartition(const Partition& partition_) -> void
    {
        // Initialise the partition member
        partition = partition_;

        // Set the partition of the equilibrium solver
        equilibrium.setPartition(partition);

        // Set the indices of the equilibrium and kinetic species
        ies = partition.indicesEquilibriumSpecies();
        iks = partition.indicesKineticSpecies();

        // Set the indices of the equilibrium and kinetic elements
        iee = partition.indicesEquilibriumElements();
        ike = partition.indicesKineticElements();

        // Set the number of equilibrium and kinetic species
        Ne = ies.size();
        Nk = iks.size();

        // Set the number of equilibrium and kinetic elements
        Ee = iee.size();
        Ek = ike.size();

        // Initialise the formula matrix of the equilibrium partition
        We = partition.formulaMatrixEquilibriumPartition();

        // Initialise the stoichiometric matrices w.r.t. the equilibrium and kinetic species
        Se = cols(reactions.stoichiometricMatrix(), ies);
        Sk = cols(reactions.stoichiometricMatrix(), iks);

        // Initialise the coefficient matrix `A` of the chemical kinetics problem
        A.resize(Ee + Nk, reactions.numReactions());
        A.topRows(Ee) = We * tr(Se);
        A.bottomRows(Nk) = tr(Sk);

        // Allocate memory for the partial derivatives of the reaction rates `r` w.r.t. to `u = [be nk]`
        drkdu.resize(reactions.numReactions(), Ee + Nk);
    }

    auto initialize(ChemicalState& state, double tstart) -> void
    {
        // Initialise the temperature and pressure variables
        T = state.temperature();
        P = state.pressure();

        // Extract the composition of the equilibrium and kinetic species
        const Vector& n = state.speciesAmounts();
        rows(n, ies).to(ne);
        rows(n, iks).to(nk);

        // Assemble the vector benk = [be nk]
        benk.resize(Ee + Nk);
        benk.segment(00, Ee) = We * ne;
        benk.segment(Ee, Nk) = nk;

        // Define the ODE function
        ODEFunction ode_function = [&](double t, const Vector& u, Vector& res)
        {
            return function(state, t, u, res);
        };

        // Define the jacobian of the ODE function
        ODEJacobian ode_jacobian = [&](double t, const Vector& u, Matrix& res)
        {
            return jacobian(state, t, u, res);
        };

        // Initialise the ODE problem
        ODEProblem problem;
        problem.setNumEquations(Ee + Nk);
        problem.setFunction(ode_function);
        problem.setJacobian(ode_jacobian);

        // Define the options for the ODE solver
        ODEOptions options_ode = options.ode;

        // Ensure the absolute tolerance values for [be, nk] take into account their initial state
        // This is important to set adequate absolute tolerances for very small components
        if(options_ode.abstols.size() == 0)
            options_ode.abstols = options_ode.abstol * benk;

        // Ensure there is no zero absolute tolerance (this causes a CVODE error)
        options_ode.abstols = (options_ode.abstols.array() <= 0).select(options_ode.abstol, options_ode.abstols);

        // Set the ODE problem and initialize the ODE solver
        ode.setProblem(problem);
        ode.setOptions(options_ode);
        ode.initialize(tstart, benk);

        // Set the options of the equilibrium solver
        equilibrium.setOptions(options.equilibrium);
    }

    auto step(ChemicalState& state, double& t) -> void
    {
        const double tfinal = unsigned(-1);
        step(state, t, tfinal);
    }

    auto step(ChemicalState& state, double& t, double tfinal) -> void
    {
        // Extract the composition vector of the equilibrium and kinetic species
        const Vector& n = state.speciesAmounts();
        rows(n, ies).to(ne);
        rows(n, iks).to(nk);

        // Assemble the vector benk = [be nk]
        benk.segment(00, Ee) = We * ne;
        benk.segment(Ee, Nk) = nk;

        // Perform one ODE step integration
        ode.integrate(t, benk, tfinal);

        // Extract the `be` and `nk` entries of the vector `benk`
        be = benk.segment(00, Ee);
        nk = benk.segment(Ee, Nk);

        // Update the composition of the kinetic species
        state.setSpeciesAmounts(nk, iks);

        // Update the composition of the equilibrium species
        equilibrium.solve(state, T, P, be);
    }

    auto solve(ChemicalState& state, double t, double dt) -> void
    {
        // Initialise the chemical kinetics solver
        initialize(state, t);

        // Integrate the chemical kinetics ODE from `t` to `t + dt`
        ode.solve(t, dt, benk);

        // Extract the `be` and `nk` entries of the vector `benk`
        be = benk.segment(00, Ee);
        nk = benk.segment(Ee, Nk);

        // Update the composition of the kinetic species
        state.setSpeciesAmounts(nk, iks);

        // Update the composition of the equilibrium species
        equilibrium.solve(state, T, P, be);
    }

    auto function(ChemicalState& state, double t, const Vector& u, Vector& res) -> int
    {
        // Extract the `be` and `nk` entries of the vector [be, nk]
        be = u.segment(00, Ee);
        nk = u.segment(Ee, Nk);

        // Check for non-finite values in the vector `benk`
        for(unsigned i = 0; i < u.rows(); ++i)
            if(!std::isfinite(u[i]))
                return 1; // ensure the ode solver will reduce the time step

        // Update the composition of the kinetic species in the member `state`
        state.setSpeciesAmounts(nk, iks);

        // Solve the equilibrium problem using the elemental molar abundance `be`
        equilibrium.solve(state, T, P, be);

        // Get the molar amounts of the species
        const Vector& n = state.speciesAmounts();

        // Update the chemical properties of the system
        properties = system.properties(T, P, n);

        // Calculate the kinetic rates of the reactions
        rk = reactions.rates(properties);

        // Calculate the right-hand side function of the ODE
        res.segment(00, Ee) = We * tr(Se) * rk.val;
        res.segment(Ee, Nk) = tr(Sk) * rk.val;

        return 0;
    }

    auto jacobian(ChemicalState& state, double t, const Vector& u, Matrix& res) -> int
    {
        // Extract the `be` and `nk` entries of the vector `benk = [be, nk]`
        be = u.segment(00, Ee);
        nk = u.segment(Ee, Nk);

        // Update the composition of the kinetic species in the member `state`
        state.setSpeciesAmounts(nk, iks);

        // Solve the equilibrium problem using the elemental molar abundance `be`
        equilibrium.solve(state, T, P, be);

        // Calculate the sensitivity of the equilibrium state
        sensitivity = equilibrium.sensitivity();

        // Extract the columns of the jacobian matrix w.r.t. the equilibrium and kinetic species
        drkdne = cols(rk.ddn, ies);
        drkdnk = cols(rk.ddn, iks);

        // Calculate the derivatives of `rk` w.r.t. `be` using the equilibrium sensitivity
        drkdbe = drkdne * sensitivity.dnedbe;

        // Assemble the partial derivatives of the reaction rates `r` w.r.t. to `u = [be nk]`
        drkdu << drkdbe, drkdnk;

        // Calculate the Jacobian matrix of the ODE function
        res = A * drkdu;

        return 0;
    }
};

KineticSolver::KineticSolver()
: pimpl(new Impl())
{}

KineticSolver::KineticSolver(const ReactionSystem& reactions)
: pimpl(new Impl(reactions))
{}

KineticSolver::KineticSolver(const KineticSolver& other)
: pimpl(new Impl(*other.pimpl))
{}

KineticSolver::~KineticSolver()
{}

auto KineticSolver::operator=(KineticSolver other) -> KineticSolver&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto KineticSolver::setOptions(const KineticOptions& options) -> void
{
    pimpl->setOptions(options);
}

auto KineticSolver::setPartition(const Partition& partition) -> void
{
    pimpl->setPartition(partition);
}

auto KineticSolver::initialize(ChemicalState& state, double tstart) -> void
{
    pimpl->initialize(state, tstart);
}

auto KineticSolver::step(ChemicalState& state, double& t) -> void
{
    pimpl->step(state, t);
}

auto KineticSolver::step(ChemicalState& state, double& t, double dt) -> void
{
    pimpl->step(state, t, dt);
}

auto KineticSolver::solve(ChemicalState& state, double t, double dt) -> void
{
    pimpl->solve(state, t, dt);
}

} // namespace Reaktoro
