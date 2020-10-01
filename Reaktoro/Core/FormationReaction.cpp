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

#include "FormationReaction.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Thermodynamics/Reactions/ReactionThermoModelConstLgK.hpp>

namespace Reaktoro {

struct FormationReaction::Impl
{
    /// The name of the product species in the formation reaction.
    String product;

    /// The reactant species in the formation reaction.
    Pairs<Species, double> reactants;

    /// The function that computes the standard thermodynamic properties of this reaction.
    ReactionThermoPropsFn rxnpropsfn;

    /// Construct a default FormationReaction::Impl object
    Impl()
    {}

    /// Return the standard thermodynamic model function of the product species.
    auto standardThermoPropsFn() const -> StandardThermoPropsFn
    {
        if(rxnpropsfn == nullptr) return {};

        return [=](real T, real P) -> StandardThermoProps
        {
            const auto R = universalGasConstant;
            StandardThermoProps props;
            ReactionThermoProps rxnprops = rxnpropsfn(T, P);

            props.G0 = rxnprops.dG0; // G0 = dG0 + sum(vr * G0r)
            props.H0 = rxnprops.dH0; // H0 = dH0 + sum(vr * H0r)
            for(const auto [reactant, coeff] : reactants) {
                const auto reactantprops = reactant.props(T, P);
                props.G0 += coeff * reactantprops.G0;
                props.H0 += coeff * reactantprops.H0;
            }
            return props;
        };
    }
};

FormationReaction::FormationReaction()
: pimpl(new Impl())
{}

auto FormationReaction::clone() const -> FormationReaction
{
    FormationReaction copy;
    *copy.pimpl = *pimpl;
    return copy;
}

auto FormationReaction::withProduct(String product) const -> FormationReaction
{
    FormationReaction copy = clone();
    copy.pimpl->product = product;
    return copy;
}

auto FormationReaction::withReactants(Pairs<Species, double> reactants) const -> FormationReaction
{
    FormationReaction copy = clone();
    copy.pimpl->reactants = reactants;
    return copy;
}

auto FormationReaction::withEquilibriumConstant(real lgK0) const -> FormationReaction
{
    return with(ReactionThermoModelConstLgK(lgK0));
}

auto FormationReaction::withReactionThermoPropsFn(const ReactionThermoPropsFn& fn) const -> FormationReaction
{
    FormationReaction copy = clone();
    copy.pimpl->rxnpropsfn = fn;
    return copy;
}

auto FormationReaction::with(const ReactionThermoPropsFn& fn) const -> FormationReaction
{
    return withReactionThermoPropsFn(fn);
}

auto FormationReaction::product() const -> String
{
    return pimpl->product;
}

auto FormationReaction::reactants() const -> const Pairs<Species, double>&
{
    return pimpl->reactants;
}

auto FormationReaction::reactionThermoPropsFn() const -> const ReactionThermoPropsFn&
{
    return pimpl->rxnpropsfn;
}

auto FormationReaction::standardThermoPropsFn() const -> StandardThermoPropsFn
{
    return pimpl->standardThermoPropsFn();
}

auto FormationReaction::stoichiometry(String reactant) const -> double
{
    for(const auto& [species, coeff] : reactants())
        if(reactant == species.name())
            return coeff;
    return 0.0;
}

} // namespace Reaktoro
