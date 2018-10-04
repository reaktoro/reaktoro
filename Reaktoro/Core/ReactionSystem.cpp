// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

#include "ReactionSystem.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/Reaction.hpp>

namespace Reaktoro {
namespace {

auto stoichiometricMatrix(const ChemicalSystem& system, const std::vector<Reaction>& reactions) -> Matrix
{
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
    /// The chemical system instance
    ChemicalSystem system;

    /// The reactions that compose the reaction system
    std::vector<Reaction> reactions;

    /// The stoichiometric matrix of the reactions w.r.t. to all species in the system
    Matrix stoichiometric_matrix;

    /// Construct a defaut ReactionSystem::Impl instance
    Impl()
    {}

    /// Construct a ReactionSystem::Impl instance with given reactions
    Impl(const ChemicalSystem& system, const std::vector<Reaction>& reactions)
    : system(system), reactions(reactions)
    {
        // Initialize the stoichiometric matrix of the reactions
        stoichiometric_matrix = Reaktoro::stoichiometricMatrix(system, reactions);
    }
};

ReactionSystem::ReactionSystem()
: pimpl(new Impl())
{}

ReactionSystem::ReactionSystem(const ChemicalSystem& system, const std::vector<Reaction>& reactions)
: pimpl(new Impl(system, reactions))
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

auto ReactionSystem::indexReactionWithError(std::string name) const -> Index
{
    const Index index = indexReaction(name);
    Assert(index < numReactions(),
        "Cannot get the index of the reaction `" + name + "`.",
        "There is no reaction called `" + name + "` in the reaction system.");
    return index;
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

auto ReactionSystem::stoichiometricMatrix() const -> MatrixConstRef
{
    return pimpl->stoichiometric_matrix;
}

auto ReactionSystem::system() const -> const ChemicalSystem&
{
    return pimpl->system;
}

auto ReactionSystem::lnEquilibriumConstants(const ChemicalProperties& properties) const -> ThermoVector
{
    const unsigned num_reactions = numReactions();
    ThermoVector res(num_reactions);
    for(unsigned i = 0; i < num_reactions; ++i)
        res[i] = reaction(i).lnEquilibriumConstant(properties);
    return res;
}

auto ReactionSystem::lnReactionQuotients(const ChemicalProperties& properties) const -> ChemicalVector
{
    const unsigned num_reactions = numReactions();
    const unsigned num_species = system().numSpecies();
    ChemicalVector res(num_reactions, num_species);
    for(unsigned i = 0; i < num_reactions; ++i)
        res[i] = reaction(i).lnReactionQuotient(properties);
    return res;
}

auto ReactionSystem::rates(const ChemicalProperties& properties) const -> ChemicalVector
{
    const unsigned num_reactions = numReactions();
    const unsigned num_species = system().numSpecies();
    ChemicalVector res(num_reactions, num_species);
    for(unsigned i = 0; i < num_reactions; ++i)
        res[i] = reaction(i).rate(properties);
    return res;
}

} // namespace Reaktoro
