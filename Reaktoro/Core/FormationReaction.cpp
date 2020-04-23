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

namespace Reaktoro {

struct FormationReaction::Impl
{
    /// The name of the product species in the formation reaction.
    String product;

    /// The reactant species in the formation reaction.
    Pairs<Species, double> reactants;

    /// The equilibrium constant function (log base 10).
    Fn<real,real,real> lgK;

    /// The enthalpy of formation function (in J/mol).
    Fn<real,real,real> dH0;

    /// Construct a default FormationReaction::Impl object
    Impl()
    {}

    /// Return the standard Gibbs energy function of the product species in the formation reaction.
    auto createStandardGibbsEnergyFn() const -> Fn<real,real,real>
    {
        return [=](real T, real P)
        {
            const auto R = universalGasConstant;
            const auto lnK = ln10 * lgK(T, P);
            const auto dG0 = -R*T*lnK;
            real G0 = dG0; // G0 = dG0 + sum(vr * G0r)
            for(const auto [reactant, coeff] : reactants) {
                assert(reactant.standardThermoPropsFn());
                G0 += coeff * reactant.standardThermoPropsFn()(T, P).G0;
            }
            return G0;
        };
    }

    /// Return the standard enthalpy function of the product species in a formation reaction.
    auto createStandardEnthalpyFn() const -> Fn<real,real,real>
    {
        return [=](real T, real P)
        {
            real H0 = dH0(T, P); // H0 = dH0 + sum(vr * H0r)
            for(const auto [reactant, coeff] : reactants) {
                assert(reactant.standardThermoPropsFn());
                H0 += coeff * reactant.standardThermoPropsFn()(T, P).H0;
            }
            return H0;
        };
    }
};

FormationReaction::FormationReaction()
: pimpl(new Impl())
{}

auto FormationReaction::withProduct(String product) const -> FormationReaction
{
    FormationReaction copy = clone();
    copy.pimpl->product = product;
    return copy;
}

auto FormationReaction::withReactants(Pairs<Species, double> reactants) const -> FormationReaction
{
    for(const auto& [species, coeff] : reactants)
        error(!species.standardThermoPropsFn(),
            "The Species objects in the method FormationReaction::withReactants "
            "need to have a non-empty function for ther standard thermodynamic "
            "property calculations. Use method Species::withStandardThermoPropsFn "
            "in the Species object with name ", species.name(), " to correct this error.");
    FormationReaction copy = clone();
    copy.pimpl->reactants = reactants;
    return copy;
}

auto FormationReaction::withEquilibriumConstant(real value) const -> FormationReaction
{
    FormationReaction copy = clone();
    copy.pimpl->lgK = [=](real T, real P) { return value; };
    return copy;
}

auto FormationReaction::withEquilibriumConstantFn(const Fn<real,real,real>& fn) const -> FormationReaction
{
    FormationReaction copy = clone();
    copy.pimpl->lgK = fn;
    return copy;
}

auto FormationReaction::withFormationEnthalpy(real value) const -> FormationReaction
{
    FormationReaction copy = clone();
    copy.pimpl->dH0 = [=](real T, real P) { return value; };
    return copy;
}

auto FormationReaction::withFormationEnthalpyFn(const Fn<real,real,real>& fn) const -> FormationReaction
{
    FormationReaction copy = clone();
    copy.pimpl->dH0 = fn;
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

auto FormationReaction::equilibriumConstantFn() const -> const Fn<real,real,real>&
{
    return pimpl->lgK;
}

auto FormationReaction::formationEnthalpyFn() const -> const Fn<real,real,real>&
{
    return pimpl->dH0;
}

auto FormationReaction::standardGibbsEnergyFn() const -> Fn<real,real,real>
{
    return pimpl->createStandardGibbsEnergyFn();
}

auto FormationReaction::standardEnthalpyFn() const -> Fn<real,real,real>
{
    return pimpl->createStandardEnthalpyFn();
}

auto FormationReaction::clone() const -> FormationReaction
{
    FormationReaction copy;
    *copy.pimpl = *pimpl;
    return copy;
}

} // namespace Reaktoro
