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

#include "KineticsSolver.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConditions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumDims.hpp>
#include <Reaktoro/Equilibrium/EquilibriumRestrictions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>
#include <Reaktoro/Kinetics/KineticsOptions.hpp>
#include <Reaktoro/Kinetics/KineticsResult.hpp>
#include <Reaktoro/Kinetics/KineticsSensitivity.hpp>
#include <Reaktoro/Kinetics/KineticsUtils.hpp>

namespace Reaktoro {

struct KineticsSolver::Impl
{
    const ChemicalSystem system;       ///< The chemical system associated with this kinetic solver.
    const EquilibriumSpecs especs;     ///< The original chemical equilibrium specifications provided at construction time.
    const EquilibriumDims edims;       ///< The original dimensions of the variables and constraints in the equilibrium specifications at construction time.
    const EquilibriumSpecs kspecs;     ///< The chemical equilibrium specifications associated with this kinetic solver.
    const EquilibriumDims kdims;       ///< The dimensions of the variables and constraints in the equilibrium specifications.
    const Index idt;                   ///< The index of the *w* input variable corresponding to Δt.
    EquilibriumSolver ksolver;         ///< The equilibrium solver used for the kinetics calculations.
    EquilibriumConditions kconditions; ///< The equilibrium conditions used for the kinetics calculations.
    KineticsOptions koptions;          ///< The options of this kinetics solver.
    KineticsResult kresult;            ///< The result of the equilibrium calculation with kinetic constraints
    VectorXr w;                        ///< The auxiliary vector used to set the input variables w of the equilibrium conditions used for the kinetics calculations.
    VectorXd c0;                       ///< The auxiliary vector used to set the initial amounts c0 of the conservative components of the equilibrium conditions used for the kinetics calculations.
    VectorXd plower;                   ///< The auxiliary vector used to set the lower bounds of p variables of the equilibrium conditions used for the kinetics calculations.
    VectorXd pupper;                   ///< The auxiliary vector used to set the upper bounds of p variables of the equilibrium conditions used for the kinetics calculations.

    /// Construct a KineticsSolver::Impl object with given equilibrium specifications to be attained during chemical kinetics.
    Impl(EquilibriumSpecs const& especs)
    : system(especs.system()),
      especs(especs),
      edims(especs),
      kspecs(detail::createEquilibriumSpecsForKinetics(especs)),
      kdims(kspecs),
      idt(kspecs.indexInputVariable("dt")),
      ksolver(kspecs),
      kconditions(kspecs),
      w(kdims.Nw),
      c0(kdims.Nc),
      plower(kdims.Np),
      pupper(kdims.Np)
    {
        // Initialize the equilibrium solver with the default options
        setOptions(koptions);
    }

    /// Set the options of the kinetics solver.
    auto setOptions(KineticsOptions const& opts) -> void
    {
        // Update the options of this kinetics solver
        koptions = opts;

        // Update the options in the underlying equilibrium solver
        ksolver.setOptions(koptions);
    }

    /// Update the equilibrium conditions for kinetics with given state and time step.
    auto updateEquilibriumConditionsForKinetics(ChemicalState& state, real const& dt) -> void
    {
        kconditions.temperature(state.temperature());
        kconditions.pressure(state.pressure());
        kconditions.setInputVariable(idt, dt);
    }

    /// Perform a kinetics step with a short time step if `state` has not reacted previously.
    auto preconditionOnFirstStep(ChemicalState& state, real const& dt) -> KineticsResult
    {
        if(state.equilibrium().empty())
        {
            updateEquilibriumConditionsForKinetics(state, koptions.dt0);
            return ksolver.solve(state, kconditions);
        }
        return {};
    }

    /// Update the equilibrium conditions for kinetics with given state, time step, and equilibrium conditions to be attained during chemical kinetics.
    auto updateEquilibriumConditionsForKinetics(ChemicalState& state, real const& dt, EquilibriumConditions const& econditions) -> void
    {
        auto const& K = system.stoichiometricMatrix();
        auto const& n0 = state.speciesAmounts();

        w << econditions.inputValues(), dt;
        c0 << econditions.initialComponentAmountsGetOrCompute(state), K.transpose() * n0.matrix();

        plower.head(edims.Np) = econditions.lowerBoundsControlVariablesP();
        plower.tail(kdims.Nr).fill(-inf); // no lower bounds for Δξ

        pupper.head(edims.Np) = econditions.upperBoundsControlVariablesP();
        pupper.tail(kdims.Nr).fill(+inf); // no upper bounds for Δξ

        kconditions.setInputVariables(w);
        kconditions.setInitialComponentAmounts(c0);
        kconditions.setLowerBoundsControlVariablesP(plower);
        kconditions.setUpperBoundsControlVariablesP(pupper);
    }

    /// Perform a kinetics step with a short time step if `state` has not reacted previously.
    auto preconditionOnFirstStep(ChemicalState& state, real const& dt, EquilibriumConditions const& econditions) -> KineticsResult
    {
        if(state.equilibrium().empty())
        {
            updateEquilibriumConditionsForKinetics(state, koptions.dt0, econditions);
            return ksolver.solve(state, kconditions);
        }
        return {};
    }

    //=================================================================================================================
    //
    // CHEMICAL KINETICS PRECONDITION METHODS
    //
    //=================================================================================================================
    auto precondition(ChemicalState& state) -> KineticsResult
    {
        return solve(state, 0.0);
    }
    auto precondition(ChemicalState& state, EquilibriumRestrictions const& restrictions) -> KineticsResult
    {
        return solve(state, 0.0, restrictions);
    }
    auto precondition(ChemicalState& state, EquilibriumConditions const& conditions) -> KineticsResult
    {
        return solve(state, 0.0, conditions);
    }
    auto precondition(ChemicalState& state, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions) -> KineticsResult
    {
        return solve(state, 0.0, conditions, restrictions);
    }

    //=================================================================================================================
    //
    // CHEMICAL KINETICS SOLVE METHODS
    //
    //=================================================================================================================

    auto solve(ChemicalState& state, real const& dt) -> KineticsResult
    {
        auto result = preconditionOnFirstStep(state, dt);
        updateEquilibriumConditionsForKinetics(state, dt);
        return result += ksolver.solve(state, kconditions);
    }

    auto solve(ChemicalState& state, real const& dt, EquilibriumRestrictions const& restrictions) -> KineticsResult
    {
        auto result = preconditionOnFirstStep(state, dt);
        updateEquilibriumConditionsForKinetics(state, dt);
        return result += ksolver.solve(state, kconditions, restrictions);
    }

    auto solve(ChemicalState& state, real const& dt, EquilibriumConditions const& conditions) -> KineticsResult
    {
        auto result = preconditionOnFirstStep(state, dt, conditions);
        updateEquilibriumConditionsForKinetics(state, dt, conditions);
        return result += ksolver.solve(state, kconditions);
    }

    auto solve(ChemicalState& state, real const& dt, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions) -> KineticsResult
    {
        auto result = preconditionOnFirstStep(state, dt, conditions);
        updateEquilibriumConditionsForKinetics(state, dt, conditions);
        return result += ksolver.solve(state, kconditions, restrictions);
    }

    //=================================================================================================================
    //
    // CHEMICAL KINETICS SOLVE METHODS WITH SENSITIVITY CALCULATION
    //
    //=================================================================================================================

    auto solve(ChemicalState& state, KineticsSensitivity& sensitivity, real const& dt) -> KineticsResult
    {
        auto result = preconditionOnFirstStep(state, dt);
        updateEquilibriumConditionsForKinetics(state, dt);
        return result += ksolver.solve(state, sensitivity, kconditions);
    }

    auto solve(ChemicalState& state, KineticsSensitivity& sensitivity, real const& dt, EquilibriumRestrictions const& restrictions) -> KineticsResult
    {
        auto result = preconditionOnFirstStep(state, dt);
        updateEquilibriumConditionsForKinetics(state, dt);
        return result += ksolver.solve(state, sensitivity, kconditions, restrictions);
    }

    auto solve(ChemicalState& state, KineticsSensitivity& sensitivity, real const& dt, EquilibriumConditions const& conditions) -> KineticsResult
    {
        auto result = preconditionOnFirstStep(state, dt, conditions);
        updateEquilibriumConditionsForKinetics(state, dt, conditions);
        return result += ksolver.solve(state, sensitivity, kconditions);
    }

    auto solve(ChemicalState& state, KineticsSensitivity& sensitivity, real const& dt, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions) -> KineticsResult
    {
        auto result = preconditionOnFirstStep(state, dt, conditions);
        updateEquilibriumConditionsForKinetics(state, dt, conditions);
        return result += ksolver.solve(state, sensitivity, kconditions, restrictions);
    }
};

KineticsSolver::KineticsSolver(ChemicalSystem const& system)
: pimpl(new Impl(EquilibriumSpecs::TP(system)))
{}

KineticsSolver::KineticsSolver(EquilibriumSpecs const& specs)
: pimpl(new Impl(specs))
{}

KineticsSolver::KineticsSolver(KineticsSolver const& other)
: pimpl(new Impl(*other.pimpl))
{}

KineticsSolver::~KineticsSolver()
{}

auto KineticsSolver::operator=(KineticsSolver other) -> KineticsSolver&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto KineticsSolver::precondition(ChemicalState& state) -> KineticsResult
{
    return pimpl->precondition(state);
}

auto KineticsSolver::precondition(ChemicalState& state, EquilibriumRestrictions const& restrictions) -> KineticsResult
{
    return pimpl->precondition(state, restrictions);
}

auto KineticsSolver::precondition(ChemicalState& state, EquilibriumConditions const& conditions) -> KineticsResult
{
    return pimpl->precondition(state, conditions);
}

auto KineticsSolver::precondition(ChemicalState& state, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions) -> KineticsResult
{
    return pimpl->precondition(state, conditions, restrictions);
}

auto KineticsSolver::solve(ChemicalState& state, real const& dt) -> KineticsResult
{
    return pimpl->solve(state, dt);
}

auto KineticsSolver::solve(ChemicalState& state, real const& dt, EquilibriumRestrictions const& restrictions) -> KineticsResult
{
    return pimpl->solve(state, dt, restrictions);
}

auto KineticsSolver::solve(ChemicalState& state, real const& dt, EquilibriumConditions const& conditions) -> KineticsResult
{
    return pimpl->solve(state, dt, conditions);
}

auto KineticsSolver::solve(ChemicalState& state, real const& dt, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions) -> KineticsResult
{
    return pimpl->solve(state, dt, conditions, restrictions);
}

auto KineticsSolver::solve(ChemicalState& state, KineticsSensitivity& sensitivity, real const& dt) -> KineticsResult
{
    return pimpl->solve(state, sensitivity, dt);
}

auto KineticsSolver::solve(ChemicalState& state, KineticsSensitivity& sensitivity, real const& dt, EquilibriumRestrictions const& restrictions) -> KineticsResult
{
    return pimpl->solve(state, sensitivity, dt, restrictions);
}

auto KineticsSolver::solve(ChemicalState& state, KineticsSensitivity& sensitivity, real const& dt, EquilibriumConditions const& conditions) -> KineticsResult
{
    return pimpl->solve(state, sensitivity, dt, conditions);
}

auto KineticsSolver::solve(ChemicalState& state, KineticsSensitivity& sensitivity, real const& dt, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions) -> KineticsResult
{
    return pimpl->solve(state, sensitivity, dt, conditions, restrictions);
}

auto KineticsSolver::setOptions(KineticsOptions const& options) -> void
{
    pimpl->setOptions(options);
}

} // namespace Reaktoro
