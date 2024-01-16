// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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

#include "KineticsUtils.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>

namespace Reaktoro {
namespace detail {

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
    // NOTE 1! There is no need to add ξ0 as inputs into `specs` for each
    // kinetically controlled reaction because these go in the b vector
    // together with element amounts.
    // NOTE 2! No need to also add surface areas as inputs because there
    // are added by default in EquilibriumSpecs implementation.
    //-------------------------------------------------------------------

    // Add Δt as input to the calculation (idt is the index of dt := Δt input in the w argument vector when defining equation constraints)
    const auto idt = specs.addInput("dt");

    // Compute matrix M = tr(K)*K
    const MatrixXd M = K.transpose() * K;

    // Add equation constraints to `specs` to model the kinetic rates of the reactions in the equilibrium problem
    EquationConstraints econstraints;
    econstraints.ids = rconstraints.ids;
    econstraints.fn = [=](ChemicalProps const& props, VectorXrConstRef const& p, VectorXrConstRef const& w) -> VectorXr
    {
        auto const& dt = w[idt]; // Δt can be found at the input vector w
        auto const& dxi = p.tail(Nr); // Δξ = the last Nr added entries in p
        const VectorXr r = props.reactionRates();
        return dxi - dt * M * r; // Δξ - ΔtMr = 0
    };

    specs.addConstraints(econstraints);

    return specs;
}

} // namespace detail
} // namespace Reaktoro
