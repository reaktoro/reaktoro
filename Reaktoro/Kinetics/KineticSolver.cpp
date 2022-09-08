// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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

#include "KineticSolver.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConditions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumDims.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProps.hpp>
#include <Reaktoro/Equilibrium/EquilibriumRestrictions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>
#include <Reaktoro/Kinetics/KineticSensitivity.hpp>

namespace Reaktoro {
namespace {

/// Return an EquilibriumSpecs object suitable for chemical kinetics
/// calculations using reactivity constraints in the equilibrium problem to
/// model the kinetic reactions.
auto createEquilibriumSpecsForKinetics(EquilibriumSpecs specs) -> EquilibriumSpecs
{
    auto const& system = specs.system();
    auto const& reactions = system.reactions();
    auto const& Nr = reactions.size();
    auto const& K = system.stoichiometricMatrix();

    // Add a p control variable to `specs` for each *extent of reaction change* variable Δξ (one for each reactive constraint, aka, restricted reaction)
    for(auto const& reaction : reactions)
    {
        ControlVariableP pvar;
        pvar.name = "dxi[" + reaction.name() + "]"; // e.g., Δξ[H2O = H+ + OH-], Δξ[Calcite]
        specs.addControlVariableP(pvar);
    }

    // The new number of p control variables in `specs` with the newly added extent of reaction variables Δξ
    auto const& Np = specs.numControlVariablesP();

    // Add reactivity constraints to `specs` to model the restricted reactions in the equilibrium problem (those that cannot react to full equilibrium but are constrained to a certain extent of reaction progress)
    ReactivityConstraints rconstraints;
    rconstraints.ids = vectorize(reactions, RKT_LAMBDA(x, x.name()));
    rconstraints.Kn = K.transpose();
    rconstraints.Kp = zeros(Nr, Np); // let p' be the last Nr entries in p; thus p' = Δξ and because tr(K)*n - Δξ = ξ0, we set Kn = tr(K) and Kp = = [0, -I] where I has dims Nr x Nr.
    rconstraints.Kp.rightCols(Nr) = -identity(Nr, Nr);

    specs.addReactivityConstraints(rconstraints);

    //-------------------------------------------------------------------
    // NOTE! There is no need to add ξ0 as inputs into `specs` for each
    // kinetically controlled reaction because these go in the b vector
    // together with element amounts
    //-------------------------------------------------------------------

    // Add surface area inputs into `specs` for each reacting phase interfaces in the system due to kinetically controlled reactions
    for(auto const& [iphase1, iphase2] : system.reactingPhaseInterfaces())
        if(iphase1 != iphase2)
            specs.addInput("SA[" + system.phase(iphase1).name() + ":" + system.phase(iphase2).name() + "]"); // SA[AqueousPhase:GaseousPhase]
        else specs.addInput("SA[" + system.phase(iphase1).name() + "]"); // SA[Calcite]

    // Add Δt as input to the calculation (idt is the index of dt := Δt input in the w argument vector when defining equation constraints)
    const auto idt = specs.addInput("dt");

    // Compute matrix M = tr(K)*K
    const MatrixXd M = K.transpose() * K;

    // Add equation constraints to `specs` to model the kinetic rates of the reactions in the equilibrium problem
    EquationConstraints econstraints;
    econstraints.ids = rconstraints.ids;
    econstraints.fn = [=](ChemicalState const& state, VectorXrConstRef p, VectorXrConstRef w) -> VectorXr
    {
        auto const& dt = w[idt]; // Δt can be found at the input vector w
        auto const& dxi = p.tail(Nr); // Δξ = the last Nr added entries in p
        const VectorXr r = state.props().reactionRates();
        return dxi - dt * M * r; // Δξ - ΔtMr = 0
    };

    specs.addConstraints(econstraints);

    return specs;
}

} // namespace

struct KineticSolver::Impl
{
    const ChemicalSystem system;       ///< The chemical system associated with this kinetic solver.
    const EquilibriumSpecs especs;     ///< The original chemical equilibrium specifications provided at construction time.
    const EquilibriumDims edims;       ///< The original dimensions of the variables and constraints in the equilibrium specifications at construction time.
    const EquilibriumSpecs kspecs;     ///< The chemical equilibrium specifications associated with this kinetic solver.
    const EquilibriumDims kdims;       ///< The dimensions of the variables and constraints in the equilibrium specifications.
    EquilibriumSolver ksolver;         ///< The equilibrium solver used for the kinetics calculations.
    EquilibriumConditions kconditions; ///< The equilibrium conditions used for the kinetics calculations.
    EquilibriumOptions options;        ///< The options of the equilibrium solver.
    EquilibriumResult result;          ///< The result of the equilibrium calculation with kinetic constraints
    VectorXr w;                        ///< The auxiliary vector used to set the inputs of the equilibrium conditions used for the kinetics calculations.
    VectorXd plower;                   ///< The auxiliary vector used to set the lower bounds of p variables of the equilibrium conditions used for the kinetics calculations.
    VectorXd pupper;                   ///< The auxiliary vector used to set the upper bounds of p variables of the equilibrium conditions used for the kinetics calculations.

    /// Construct a Impl instance with given EquilibriumConditions object.
    Impl(EquilibriumSpecs const& especs)
    : system(especs.system()),
      especs(especs),
      edims(especs),
      kspecs(createEquilibriumSpecsForKinetics(especs)),
      kdims(kspecs),
      ksolver(kspecs),
      kconditions(kspecs),
      w(kdims.Nw),
      plower(kdims.Np),
      pupper(kdims.Np)
    {
        // Initialize the equilibrium solver with the default options
        setOptions(options);
    }

    /// Set the options of the equilibrium solver.
    auto setOptions(EquilibriumOptions const& opts) -> void
    {
        // Update the options of the equilibrium calculation
        options = opts;

        // Update the options in the underlying equilibrium solver
        ksolver.setOptions(options);

        // // Ensure some options have proper values
        // error(options.epsilon <= 0, "EquilibriumOptions::epsilon cannot be zero or negative.");

        // // Initialize the names of the primal and dual variables
        // if(options.optima.output.active)
        // {
        //     // Define some auxiliary references to the variables names
        //     auto& xnames = options.optima.output.xnames;

        //     // Initialize the names of the variables corresponding to the species
        //     for(auto species : system.species())
        //         xnames.push_back("n[" + species.name() + "]");

        //     // Initialize the names of the variables corresponding to the implicit titrants
        //     for(auto titrant : specs.namesTitrantsImplicit())
        //         xnames.push_back(titrant);
        // }
    }

    auto updateEquilibriumConditions(ChemicalState& state, real const& dt, EquilibriumConditions const& econditions) -> void
    {
        w << econditions.inputValues(), state.surfaceAreas(), dt;

        plower.head(edims.Np) = econditions.lowerBoundsControlVariablesP();
        plower.tail(kdims.Nr).fill(-inf); // no lower bounds for Δξ

        pupper.head(edims.Np) = econditions.upperBoundsControlVariablesP();
        pupper.tail(kdims.Nr).fill(+inf); // no upper bounds for Δξ
    }

    //=================================================================================================================
    //
    // CHEMICAL KINETICS METHODS
    //
    //=================================================================================================================

    auto solve(ChemicalState& state, real const& dt) -> EquilibriumResult
    {
        EquilibriumRestrictions restrictions(system);
        return solve(state, dt, restrictions);
    }

    auto solve(ChemicalState& state, real const& dt, EquilibriumRestrictions const& restrictions) -> EquilibriumResult
    {
        EquilibriumConditions conditions(especs);
        conditions.temperature(state.temperature());
        conditions.pressure(state.pressure());
        return solve(state, dt, conditions, restrictions);
    }

    auto solve(ChemicalState& state, real const& dt, EquilibriumConditions const& conditions) -> EquilibriumResult
    {
        EquilibriumRestrictions restrictions(system);
        return solve(state, dt, conditions, restrictions);
    }

    auto solve(ChemicalState& state, real const& dt, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions) -> EquilibriumResult
    {
        updateEquilibriumConditions(state, dt, conditions);
        return ksolver.solve(state, kconditions, restrictions);
    }

    //=================================================================================================================
    //
    // CHEMICAL KINETICS METHODS WITH SENSITIVITY CALCULATION
    //
    //=================================================================================================================

    auto solve(ChemicalState& state, KineticSensitivity& sensitivity, real const& dt) -> EquilibriumResult
    {
        EquilibriumRestrictions restrictions(system);
        return solve(state, sensitivity, dt, restrictions);
    }

    auto solve(ChemicalState& state, KineticSensitivity& sensitivity, real const& dt, EquilibriumRestrictions const& restrictions) -> EquilibriumResult
    {
        EquilibriumConditions conditions(especs);
        conditions.temperature(state.temperature());
        conditions.pressure(state.pressure());
        return solve(state, sensitivity, dt, conditions, restrictions);
    }

    auto solve(ChemicalState& state, KineticSensitivity& sensitivity, real const& dt, EquilibriumConditions const& conditions) -> EquilibriumResult
    {
        EquilibriumRestrictions restrictions(system);
        return solve(state, sensitivity, dt, conditions, restrictions);
    }

    auto solve(ChemicalState& state, KineticSensitivity& sensitivity, real const& dt, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions) -> EquilibriumResult
    {
        updateEquilibriumConditions(state, dt, conditions);
        return ksolver.solve(state, sensitivity, kconditions, restrictions);
    }

    //=================================================================================================================
    //
    // CHEMICAL KINETICS METHODS WITH GIVEN AMOUNTS OF CONSERVATIVE COMPONENTS
    //
    //=================================================================================================================

    auto solve(ChemicalState& state, real const& dt, ArrayXdConstRef c0) -> EquilibriumResult
    {
        EquilibriumRestrictions restrictions(system);
        return solve(state, dt, restrictions, c0);
    }

    auto solve(ChemicalState& state, real const& dt, EquilibriumRestrictions const& restrictions, ArrayXdConstRef c0) -> EquilibriumResult
    {
        EquilibriumConditions conditions(especs);
        conditions.temperature(state.temperature());
        conditions.pressure(state.pressure());
        return solve(state, dt, conditions, restrictions, c0);
    }

    auto solve(ChemicalState& state, real const& dt, EquilibriumConditions const& conditions, ArrayXdConstRef c0) -> EquilibriumResult
    {
        EquilibriumRestrictions restrictions(system);
        return solve(state, dt, conditions, restrictions, c0);
    }

    auto solve(ChemicalState& state, real const& dt, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions, ArrayXdConstRef c0) -> EquilibriumResult
    {
        updateEquilibriumConditions(state, dt, conditions);
        return ksolver.solve(state, kconditions, restrictions, c0);
    }

    //=================================================================================================================
    //
    // CHEMICAL KINETICS METHODS WITH GIVEN AMOUNTS OF CONSERVATIVE COMPONENTS AND SENSITIVITY CALCULATION
    //
    //=================================================================================================================

    auto solve(ChemicalState& state, KineticSensitivity& sensitivity, real const& dt, ArrayXdConstRef c0) -> EquilibriumResult
    {
        EquilibriumRestrictions restrictions(system);
        return solve(state, sensitivity, dt, restrictions, c0);
    }

    auto solve(ChemicalState& state, KineticSensitivity& sensitivity, real const& dt, EquilibriumRestrictions const& restrictions, ArrayXdConstRef c0) -> EquilibriumResult
    {
        EquilibriumConditions conditions(especs);
        conditions.temperature(state.temperature());
        conditions.pressure(state.pressure());
        return solve(state, sensitivity, dt, conditions, restrictions, c0);
    }

    auto solve(ChemicalState& state, KineticSensitivity& sensitivity, real const& dt, EquilibriumConditions const& conditions, ArrayXdConstRef c0) -> EquilibriumResult
    {
        EquilibriumRestrictions restrictions(system);
        return solve(state, sensitivity, dt, conditions, restrictions, c0);
    }

    auto solve(ChemicalState& state, KineticSensitivity& sensitivity, real const& dt, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions, ArrayXdConstRef c0) -> EquilibriumResult
    {
        updateEquilibriumConditions(state, dt, conditions);
        return ksolver.solve(state, sensitivity, kconditions, restrictions, c0);
    }
};

KineticSolver::KineticSolver(ChemicalSystem const& system)
: pimpl(new Impl(EquilibriumSpecs::TP(system)))
{}

KineticSolver::KineticSolver(EquilibriumSpecs const& specs)
: pimpl(new Impl(specs))
{}

KineticSolver::KineticSolver(KineticSolver const& other)
: pimpl(new Impl(*other.pimpl))
{}

KineticSolver::~KineticSolver()
{}

auto KineticSolver::operator=(KineticSolver other) -> KineticSolver&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto KineticSolver::solve(ChemicalState& state, real const& dt) -> EquilibriumResult
{
    return pimpl->solve(state, dt);
}

auto KineticSolver::solve(ChemicalState& state, real const& dt, EquilibriumRestrictions const& restrictions) -> EquilibriumResult
{
    return pimpl->solve(state, dt, restrictions);
}

auto KineticSolver::solve(ChemicalState& state, real const& dt, EquilibriumConditions const& conditions) -> EquilibriumResult
{
    return pimpl->solve(state, dt, conditions);
}

auto KineticSolver::solve(ChemicalState& state, real const& dt, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions) -> EquilibriumResult
{
    return pimpl->solve(state, dt, conditions, restrictions);
}

auto KineticSolver::solve(ChemicalState& state, KineticSensitivity& sensitivity, real const& dt) -> EquilibriumResult
{
    return pimpl->solve(state, sensitivity, dt);
}

auto KineticSolver::solve(ChemicalState& state, KineticSensitivity& sensitivity, real const& dt, EquilibriumRestrictions const& restrictions) -> EquilibriumResult
{
    return pimpl->solve(state, sensitivity, dt, restrictions);
}

auto KineticSolver::solve(ChemicalState& state, KineticSensitivity& sensitivity, real const& dt, EquilibriumConditions const& conditions) -> EquilibriumResult
{
    return pimpl->solve(state, sensitivity, dt, conditions);
}

auto KineticSolver::solve(ChemicalState& state, KineticSensitivity& sensitivity, real const& dt, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions) -> EquilibriumResult
{
    return pimpl->solve(state, sensitivity, dt, conditions, restrictions);
}

auto KineticSolver::solve(ChemicalState& state, real const& dt, ArrayXdConstRef c0) -> EquilibriumResult
{
    return pimpl->solve(state, dt, c0);
}

auto KineticSolver::solve(ChemicalState& state, real const& dt, EquilibriumRestrictions const& restrictions, ArrayXdConstRef c0) -> EquilibriumResult
{
    return pimpl->solve(state, dt, restrictions, c0);
}

auto KineticSolver::solve(ChemicalState& state, real const& dt, EquilibriumConditions const& conditions, ArrayXdConstRef c0) -> EquilibriumResult
{
    return pimpl->solve(state, dt, conditions, c0);
}

auto KineticSolver::solve(ChemicalState& state, real const& dt, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions, ArrayXdConstRef c0) -> EquilibriumResult
{
    return pimpl->solve(state, dt, conditions, restrictions, c0);
}

auto KineticSolver::solve(ChemicalState& state, KineticSensitivity& sensitivity, real const& dt, ArrayXdConstRef c0) -> EquilibriumResult
{
    return pimpl->solve(state, sensitivity, dt, c0);
}

auto KineticSolver::solve(ChemicalState& state, KineticSensitivity& sensitivity, real const& dt, EquilibriumRestrictions const& restrictions, ArrayXdConstRef c0) -> EquilibriumResult
{
    return pimpl->solve(state, sensitivity, dt, restrictions, c0);
}

auto KineticSolver::solve(ChemicalState& state, KineticSensitivity& sensitivity, real const& dt, EquilibriumConditions const& conditions, ArrayXdConstRef c0) -> EquilibriumResult
{
    return pimpl->solve(state, sensitivity, dt, conditions, c0);
}

auto KineticSolver::solve(ChemicalState& state, KineticSensitivity& sensitivity, real const& dt, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions, ArrayXdConstRef c0) -> EquilibriumResult
{
    return pimpl->solve(state, sensitivity, dt, conditions, restrictions, c0);
}

auto KineticSolver::setOptions(EquilibriumOptions const& options) -> void
{
    pimpl->setOptions(options);
}

} // namespace Reaktoro
