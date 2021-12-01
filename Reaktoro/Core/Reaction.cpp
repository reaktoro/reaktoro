// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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

#include "Reaction.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Core/SpeciesList.hpp>
#include <Reaktoro/Core/Utils.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>

namespace Reaktoro {

struct Reaction::Impl
{
    /// The name that uniquely identifies this reaction.
    String name;

    /// The equation of the reaction with its species and stoichiometric coefficients.
    ReactionEquation equation;

    /// The chemical system instance
    ChemicalSystem system;

    /// The function that computes the equilibrium constant of the reaction (in natural log).
    EquilibriumConstantFn lnKfn;

    /// The function that computes the rate of the reaction (in mol/s).
    ReactionRateFn ratefn;

    /// The species in the reaction
    SpeciesList species;

    /// The indices of the species in the reaction
    Indices indices;

    /// The stoichiometries of the species in the reaction
    ArrayXr stoichiometries;

    Impl()
    {}

    Impl(const ReactionEquation& equation, const ChemicalSystem& system)
    : equation(equation), system(system)
    {
        auto numspecies = equation.species().size();

        // Initialize the species, their indices, and their stoichiometries in the reaction
        stoichiometries.resize(numspecies);
        unsigned i = 0;
        for(const auto& pair : equation.equation())
        {
            Index idx = detail::resolveSpeciesIndex(system, pair.first.name());
            species.push_back(pair.first);
            indices.push_back(idx);
            stoichiometries[i] = pair.second;
            ++i;
        }
    }

};

Reaction::Reaction()
: pimpl(new Impl())
{}

auto Reaction::clone() const -> Reaction
{
    Reaction reaction;
    *reaction.pimpl = *pimpl;
    return reaction;
}

auto Reaction::withName(String name) const -> Reaction
{
    Reaction copy = clone();
    copy.pimpl->name = name;
    return copy;
}

auto Reaction::withEquation(const ReactionEquation& equation) const -> Reaction
{
    Reaction copy = clone();
    copy.pimpl->equation = equation;
    return copy;
}

auto Reaction::withEquilibriumConstantFn(const EquilibriumConstantFn& fn) const -> Reaction
{
    Reaction copy = clone();
    copy.pimpl->lnKfn = fn;
    return copy;
}

auto Reaction::withRateFn(const ReactionRateFn& fn) const -> Reaction
{
    Reaction copy = clone();
    copy.pimpl->ratefn = fn;
    return copy;
}

auto Reaction::name() const -> String
{
    return pimpl->name;
}

auto Reaction::equation() const -> const ReactionEquation&
{
    return pimpl->equation;
}

auto Reaction::equilibriumConstantFn() const -> const EquilibriumConstantFn&
{
    return pimpl->lnKfn;
}

auto Reaction::rateFn() const -> const ReactionRateFn&
{
    return pimpl->ratefn;
}

auto Reaction::lnEquilibriumConstant(const ChemicalProps& properties) const -> real
{
    // Get the temperature and pressure of the system
    const auto T = properties.temperature();
    const auto P = properties.pressure();

    // Check if a equilibrium constant function was provided
    if(pimpl->lnKfn) return pimpl->lnKfn(T, P);

    // Calculate the equilibrium constant using the standard Gibbs energies of the species
    const auto G0 = properties.speciesStandardGibbsEnergies();
    const auto RT = universalGasConstant * T;

    real res;
    for(unsigned i = 0; i < pimpl->indices.size(); ++i)
    {
        const Index ispecies = pimpl->indices[i];
        const double vi = pimpl->stoichiometries[i];
        const auto G0i = G0[ispecies];
        res += vi * G0i;
    }

    return -res/RT;
}

auto Reaction::lnReactionQuotient(const ChemicalProps& properties) const -> real
{
    const auto& ln_a = properties.speciesActivitiesLn();
    real ln_Q = 0.0;
    unsigned counter = 0;
    for(Index i : pimpl->indices)
    {
        const auto vi = pimpl->stoichiometries[counter];
        ln_Q += vi * ln_a[i].val();
        ++counter;
    }

    return ln_Q;
}

auto Reaction::lnEquilibriumIndex(const ChemicalProps& properties) const -> real
{
    return lnReactionQuotient(properties) - lnEquilibriumConstant(properties);
}

auto Reaction::stoichiometry(std::string species) const -> double
{
    return pimpl->equation.stoichiometry(species);
}

auto operator<(const Reaction& lhs, const Reaction& rhs) -> bool
{
    return lhs.name() < rhs.name();
}

auto operator==(const Reaction& lhs, const Reaction& rhs) -> bool
{

    return lhs.name() == rhs.name();
}

} // namespace Reaktoro
