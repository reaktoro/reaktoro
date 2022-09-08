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
    for(auto i = 0; i < reactions.size(); ++i)
    {
        ControlVariableP pvar;
        pvar.name = "dxi[" + reactions[i].name() + "]"; // e.g., Δξ[H2O = H+ + OH-], Δξ[Calcite]
        specs.addControlVariableP(pvar);
    }

    // The new number of p control variables in `specs` with the newly added extent of reaction variables Δξ
    auto const& Np = specs.numControlVariablesP();

    // The offset that separates original p variables from the new ones
    auto const offset = Np - Nr;

    // Add reactivity constraints to `specs` to model the restricted reactions in the equilibrium problem (those that cannot react to full equilibrium but are constrained to a certain extent of reaction progress)
    for(auto i = 0; i < K.cols(); ++i)
    {
        ReactivityConstraint rconstraint;
        rconstraint.id = reactions[i].name();
        rconstraint.Kn = K.col(i);
        rconstraint.Kp = -unit(Np, offset + i); // let p' be the last Nr entries in p; thus p' = Δξ and because tr(K)*n - Δξ = ξ0, we use negative unit here so that we obtain this equation with tr(Kn)*n + tr(Kp)*p = ξ0.
        specs.addReactivityConstraint(rconstraint);
    }

    const auto idt = specs.addInput("dt"); // Add dt as input to the calculation (idt is the index of dt input in the w argument vector when defining equation constraints)

    // Add equation constraints to `specs` to model the kinetic rates of the reactions in the equilibrium problem
    for(auto i = 0; i < K.cols(); ++i)
    {
        ConstraintEquation econstraint;
        econstraint.id = reactions[i].name();
        const auto idx = specs.addInput("xi0[" + econstraint.id + "]");
        econstraint.fn = [=](ChemicalState const& state, VectorXrConstRef p, VectorXrConstRef w) -> real
        {
            auto const& dt = w[idt];
            auto const& pi = p[offset + i];
            // return pi - dt * M.row(i) * r;
            return {};
        };

        specs.addConstraint(econstraint);
    }


    return specs;
}

} // namespace

struct KineticSolver::Impl
{
    /// The chemical equilibrium specifications associated with this kinetic solver.
    const EquilibriumSpecs specs;

    /// The chemical system associated with this kinetic solver.
    const ChemicalSystem system;

    /// The dimensions of the variables and constraints in the equilibrium specifications.
    const EquilibriumDims dims;

    /// The equilibrium solver used for the kinetics calculations.
    const EquilibriumSolver solver;

    /// The options of the equilibrium solver.
    EquilibriumOptions options;

    // The equilibrium result
    EquilibriumResult result;

    /// Construct a Impl instance with given EquilibriumConditions object.
    Impl(EquilibriumSpecs const& especs)
    : specs(createEquilibriumSpecsForKinetics(especs)), system(specs.system()), dims(specs), solver(specs)
    {
        // Initialize the equilibrium solver with the default options
        setOptions(options);
    }

    /// Set the options of the equilibrium solver.
    auto setOptions(EquilibriumOptions const& opts) -> void
    {
        // Update the options of the equilibrium calculation
        options = opts;

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

    auto solve(ChemicalState& state) -> EquilibriumResult
    {
        return {};
    }

    auto solve(ChemicalState& state, EquilibriumRestrictions const& restrictions) -> EquilibriumResult
    {
        return {};
    }

    auto solve(ChemicalState& state, EquilibriumConditions const& conditions) -> EquilibriumResult
    {
        return {};
    }

    auto solve(ChemicalState& state, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions) -> EquilibriumResult
    {
        return {};
    }

    auto solve(ChemicalState& state, KineticSensitivity& sensitivity) -> EquilibriumResult
    {
        return {};
    }

    auto solve(ChemicalState& state, KineticSensitivity& sensitivity, EquilibriumRestrictions const& restrictions) -> EquilibriumResult
    {
        return {};
    }

    auto solve(ChemicalState& state, KineticSensitivity& sensitivity, EquilibriumConditions const& conditions) -> EquilibriumResult
    {
        return {};
    }

    auto solve(ChemicalState& state, KineticSensitivity& sensitivity, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions) -> EquilibriumResult
    {
        return {};
    }

    auto solve(ChemicalState& state, ArrayXdConstRef b0) -> EquilibriumResult
    {
        return {};
    }

    auto solve(ChemicalState& state, EquilibriumRestrictions const& restrictions, ArrayXdConstRef b0) -> EquilibriumResult
    {
        return {};
    }

    auto solve(ChemicalState& state, EquilibriumConditions const& conditions, ArrayXdConstRef b0) -> EquilibriumResult
    {
        return {};
    }

    auto solve(ChemicalState& state, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions, ArrayXdConstRef b0) -> EquilibriumResult
    {
        return {};
    }

    auto solve(ChemicalState& state, KineticSensitivity& sensitivity, ArrayXdConstRef b0) -> EquilibriumResult
    {
        return {};
    }

    auto solve(ChemicalState& state, KineticSensitivity& sensitivity, EquilibriumRestrictions const& restrictions, ArrayXdConstRef b0) -> EquilibriumResult
    {
        return {};
    }

    auto solve(ChemicalState& state, KineticSensitivity& sensitivity, EquilibriumConditions const& conditions, ArrayXdConstRef b0) -> EquilibriumResult
    {
        return {};
    }

    auto solve(ChemicalState& state, KineticSensitivity& sensitivity, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions, ArrayXdConstRef b0) -> EquilibriumResult
    {
        return {};
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

auto KineticSolver::solve(ChemicalState& state) -> EquilibriumResult
{
    return pimpl->solve(state);
}

auto KineticSolver::solve(ChemicalState& state, EquilibriumRestrictions const& restrictions) -> EquilibriumResult
{
    return pimpl->solve(state, restrictions);
}

auto KineticSolver::solve(ChemicalState& state, EquilibriumConditions const& conditions) -> EquilibriumResult
{
    return pimpl->solve(state, conditions);
}

auto KineticSolver::solve(ChemicalState& state, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions) -> EquilibriumResult
{
    return pimpl->solve(state, conditions, restrictions);
}

auto KineticSolver::solve(ChemicalState& state, KineticSensitivity& sensitivity) -> EquilibriumResult
{
    return pimpl->solve(state, sensitivity);
}

auto KineticSolver::solve(ChemicalState& state, KineticSensitivity& sensitivity, EquilibriumRestrictions const& restrictions) -> EquilibriumResult
{
    return pimpl->solve(state, sensitivity, restrictions);
}

auto KineticSolver::solve(ChemicalState& state, KineticSensitivity& sensitivity, EquilibriumConditions const& conditions) -> EquilibriumResult
{
    return pimpl->solve(state, sensitivity, conditions);
}

auto KineticSolver::solve(ChemicalState& state, KineticSensitivity& sensitivity, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions) -> EquilibriumResult
{
    return pimpl->solve(state, sensitivity, conditions, restrictions);
}

auto KineticSolver::solve(ChemicalState& state, ArrayXdConstRef b0) -> EquilibriumResult
{
    return pimpl->solve(state, b0);
}

auto KineticSolver::solve(ChemicalState& state, EquilibriumRestrictions const& restrictions, ArrayXdConstRef b0) -> EquilibriumResult
{
    return pimpl->solve(state, restrictions, b0);
}

auto KineticSolver::solve(ChemicalState& state, EquilibriumConditions const& conditions, ArrayXdConstRef b0) -> EquilibriumResult
{
    return pimpl->solve(state, conditions, b0);
}

auto KineticSolver::solve(ChemicalState& state, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions, ArrayXdConstRef b0) -> EquilibriumResult
{
    return pimpl->solve(state, conditions, restrictions, b0);
}

auto KineticSolver::solve(ChemicalState& state, KineticSensitivity& sensitivity, ArrayXdConstRef b0) -> EquilibriumResult
{
    return pimpl->solve(state, sensitivity, b0);
}

auto KineticSolver::solve(ChemicalState& state, KineticSensitivity& sensitivity, EquilibriumRestrictions const& restrictions, ArrayXdConstRef b0) -> EquilibriumResult
{
    return pimpl->solve(state, sensitivity, restrictions, b0);
}

auto KineticSolver::solve(ChemicalState& state, KineticSensitivity& sensitivity, EquilibriumConditions const& conditions, ArrayXdConstRef b0) -> EquilibriumResult
{
    return pimpl->solve(state, sensitivity, conditions, b0);
}

auto KineticSolver::solve(ChemicalState& state, KineticSensitivity& sensitivity, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions, ArrayXdConstRef b0) -> EquilibriumResult
{
    return pimpl->solve(state, sensitivity, conditions, restrictions, b0);
}

auto KineticSolver::setOptions(EquilibriumOptions const& options) -> void
{
    pimpl->setOptions(options);
}


} // namespace Reaktoro
