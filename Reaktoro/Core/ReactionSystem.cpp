// Reaktoro is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
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

#include "ReactionSystem.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ChemicalSystemProperties.hpp>
#include <Reaktoro/Core/Reaction.hpp>
#include <Reaktoro/Core/Utils.hpp>

namespace Reaktoro {
namespace {

auto stoichiometricMatrix(const std::vector<Reaction>& reactions) -> Matrix
{
    const auto& system = reactions.front().system();
    const auto& species = system.species();
    const auto num_reactions = reactions.size();
    const auto num_species = species.size();
    Matrix S = zeros(num_reactions, num_species);
    for(unsigned i = 0; i < num_reactions; ++i)
        for(unsigned j = 0; j < num_species; ++j)
            S(i, j) = reactions[i].stoichiometry(species[j].name());
    return S;
}

} // namespace

struct ReactionSystem::Impl
{
    /// The reactions that compose the reaction system
    std::vector<Reaction> reactions;

    /// The chemical system instance
    ChemicalSystem system;

    /// The stoichiometric matrix of the reactions w.r.t. to all species in the system
    Matrix stoichiometric_matrix;

    /// The function for the equilibrium constant of the reactions (in natural log).
    ThermoVectorFunction lnk;

    /// The function for the molar volumes of the reactions (in units of m3/mol).
    ReactionRateVectorFunction rates;

    /// Construct a defaut ReactionSystem::Impl instance
    Impl()
    {}

    /// Construct a ReactionSystem::Impl instance with given reactions
    Impl(const std::vector<Reaction>& _reactions)
    : reactions(_reactions)
    {
        // Assert the given reactions are not empty
        Assert(reactions.size(),
            "Cannot construct the ReactionSystem instance with given reactions.",
            "The given collection of reactions are empty.");

        // Initialize the systemhemical system instance
        system = reactions.front().system();

        // Initialize the stoichiometric matrix of the reactions
        stoichiometric_matrix = Reaktoro::stoichiometricMatrix(reactions);

        const unsigned num_species = system.numSpecies();
        const unsigned num_reactions = reactions.size();

        lnk = [=](double T, double P)
        {
            ThermoVector res(num_reactions);
            for(unsigned i = 0; i < num_reactions; ++i)
                res.row(i) = reactions[i].lnEquilibriumConstant(T, P);
            return res;
        };

        rates = [=](const ChemicalSystemProperties& properties)
        {
            ChemicalVector res(num_reactions, num_species);
            for(unsigned i = 0; i < num_reactions; ++i)
                res.row(i) = reactions[i].rate(properties);
            return res;
        };
    }
};

ReactionSystem::ReactionSystem()
: pimpl(new Impl())
{}

ReactionSystem::ReactionSystem(const std::vector<Reaction>& reactions)
: pimpl(new Impl(reactions))
{}

ReactionSystem::~ReactionSystem()
{}

auto ReactionSystem::numReactions() const -> unsigned
{
    return reactions().size();
}

auto ReactionSystem::indexReaction(std::string name) const -> Index
{
    return index(name, reactions());
}

auto ReactionSystem::reactions() const -> const std::vector<Reaction>&
{
    return pimpl->reactions;
}

auto ReactionSystem::reaction(Index index) const -> const Reaction&
{
    Assert(index < numReactions(),
        "Cannot return a Reaction instance with given "
        "index `" + std::to_string(index) + "`.",
        "The reaction index must be less than the "
        "number of reactions `" + std::to_string(numReactions()) + "`.");

    return pimpl->reactions[index];
}

auto ReactionSystem::reaction(std::string name) const -> const Reaction&
{
    const Index index = indexReaction(name);

    Assert(index < numReactions(),
        "Cannot return a Reaction instance with given name `" + name + "`.",
        "There is no reaction with such name in this ReactionSystem instance.");

    return pimpl->reactions[index];
}

auto ReactionSystem::stoichiometricMatrix() const -> const Matrix&
{
    return pimpl->stoichiometric_matrix;
}

auto ReactionSystem::system() const -> const ChemicalSystem&
{
    return pimpl->system;
}

auto ReactionSystem::lnEquilibriumConstants(double T, double P) const -> ThermoVector
{
    return pimpl->lnk(T, P);
}

auto ReactionSystem::rates(const ChemicalSystemProperties& properties) const -> ChemicalVector
{
    return pimpl->rates(properties);
}

auto ReactionSystem::lnReactionQuotients(const ChemicalSystemProperties& properties) const -> ChemicalVector
{
    const unsigned num_reactions = pimpl->reactions.size();
    const unsigned num_species = pimpl->system.numSpecies();
    ChemicalVector res(num_reactions, num_species);
    for(unsigned i = 0; i < num_reactions; ++i)
        res.row(i) = reaction(i).lnReactionQuotient(properties);
    return res;
}

} // namespace Reaktoro
