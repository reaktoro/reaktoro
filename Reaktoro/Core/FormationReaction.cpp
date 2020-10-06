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

    /// The function that computes the standard molar volume of the product species.
    Model<real(real,real)> std_volume_model;

    /// The function that computes the standard thermodynamic properties of this reaction.
    ReactionThermoModel rxn_thermo_model;

    /// Construct a default FormationReaction::Impl object
    Impl()
    {}

    /// Return the standard thermodynamic model function of the product species.
    auto standardThermoPropsFn() const -> StandardThermoPropsFn
    {
        error(!rxn_thermo_model.initialized(), "Could not create the standard thermodynamic "
            "model function of species ", product, " because no reaction thermodynamic "
            "model has been set. Use one of the methods below to correct this: \n"
            "    1) FormationReaction::withEquilibriumConstant\n"
            "    2) FormationReaction::withReactionThermoModel");

        error(!std_volume_model.initialized(), "Could not create the standard thermodynamic "
            "model function of species ", product, " because no standard molar volume "
            "constant or function has been set. Use one of the methods below to correct this: \n"
            "    1) FormationReaction::withProductStandardVolume\n"
            "    2) FormationReaction::withProductStandardVolumeModel");

        const auto num_reactants = reactants.size();

        Vec<StandardThermoProps> reactants_props(num_reactants);

        return [=](real T, real P) mutable -> StandardThermoProps
        {
            // Precompute the standard thermo properties of each reactant species
            for(auto i = 0; i < num_reactants; ++i)
            {
                const auto& reactant = reactants[i].first;
                reactants_props[i] = reactant.props(T, P);
            }

            // Compute the standard molar volume of the product species
            const auto V0p = std_volume_model(T, P);

            // Compute the standard molar volume change of the reaction
            auto dV0 = V0p;
            for(auto i = 0; i < num_reactants; ++i)
            {
                const auto& coeff = reactants[i].second;
                dV0 -= coeff * reactants_props[i].V0; // coeff is positve for left-hand side reactant, negative for right-hand side
            }

            // Compute the rest of the standard thermodynamic properties of the reaction
            ReactionThermoProps rxnprops;
            rxn_thermo_model.apply(rxnprops, {T, P, dV0});

            // Compute finally the standard thermodynamic properties of the product species
            StandardThermoProps props;

            props.V0 = V0p;
            props.G0 = rxnprops.dG0; // G0 = dG0 + sum(vr * G0r)
            props.H0 = rxnprops.dH0; // H0 = dH0 + sum(vr * H0r)
            for(auto i = 0; i < num_reactants; ++i)
            {
                const auto& coeff = reactants[i].second;
                const auto& reactantprops = reactants_props[i];
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
    FormationReaction copy = clone();
    copy = copy.withReactionThermoModel(ReactionThermoModelConstLgK(lgK0));
    copy = copy.withProductStandardVolume(0.0);
    return copy;
}

auto FormationReaction::withProductStandardVolume(real V0p) const -> FormationReaction
{
    return withProductStandardVolumeModel(Model<real(real,real)>::Constant(V0p, "V0"));
}

auto FormationReaction::withProductStandardVolumeModel(Model<real(real,real)> fn) const -> FormationReaction
{
    FormationReaction copy = clone();
    copy.pimpl->std_volume_model = fn;
    return copy;
}

auto FormationReaction::withReactionThermoModel(const ReactionThermoModel& fn) const -> FormationReaction
{
    FormationReaction copy = clone();
    copy.pimpl->rxn_thermo_model = fn;
    return copy;
}

auto FormationReaction::product() const -> String
{
    return pimpl->product;
}

auto FormationReaction::reactants() const -> const Pairs<Species, double>&
{
    return pimpl->reactants;
}

auto FormationReaction::reactionThermoModel() const -> const ReactionThermoModel&
{
    return pimpl->rxn_thermo_model;
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
