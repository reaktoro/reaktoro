/*
 * Reaktor is a C++ library for computational reaction modelling.
 *
 * Copyright (C) 2014 Allan Leal
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "ReactionSystem.hpp"

// C++ includes
#include <cmath>
#include <iomanip>
#include <set>

// Reaktor includes
#include <Reaktor/Core/ChemicalState.hpp>
#include <Reaktor/Core/ChemicalSystem.hpp>
#include <Reaktor/Core/Species.hpp>
#include <Reaktor/Utils/SetUtils.hpp>

namespace Reaktor {

ReactionSystem::ReactionSystem()
{}

ReactionSystem::ReactionSystem(const std::vector<Reaction>& reactions, const ChemicalSystem& system)
: m_reactions(reactions)
{
    initialiseStoichiometricMatrix(system);
    initialiseReactionPhaseMap(system);
    initialiseReactionSpeciesMap();
    initialiseSpeciesReactionMap(system);
}

auto ReactionSystem::numReactions() const -> unsigned
{
    return m_reactions.size();
}

auto ReactionSystem::reactions() const -> const std::vector<Reaction>&
{
    return m_reactions;
}

auto ReactionSystem::stoichiometricMatrix() const -> const Matrix&
{
    return stoichiometric_matrix;
}

auto ReactionSystem::mapReactionToPhase() const -> const std::vector<Indices>&
{
    return reaction_phase_map;
}

auto ReactionSystem::mapReactionToSpecies() const -> const std::vector<Indices>&
{
    return reaction_species_map;
}

auto ReactionSystem::mapSpeciesToReaction() const -> const std::vector<Indices>&
{
    return species_reaction_map;
}

auto ReactionSystem::equilibriumConstants(double T, double P) const -> Vector
{
    Vector K(m_reactions.size());

    for(int r = 0; r < K.rows(); ++r)
        K[r] = m_reactions[r].equilibriumConstant(T, P);

    return K;
}

auto ReactionSystem::reactionQuotients(const VectorResult& a) const -> VectorResult
{
    const unsigned R = m_reactions.size();
    const unsigned N = func(a).size();

    VectorResult res;
    func(res) = zeros(R);
    grad(res) = zeros(R, N);

    ScalarResult aux;
    for(unsigned i = 0; i < R; ++i)
    {
        aux = m_reactions[i].reactionQuotient(a);

        func(res)[i]     = func(aux);
        grad(res).row(i) = grad(aux);
    }

    return res;
}

auto ReactionSystem::rates(double T, double P, const Vector& n, const VectorResult& a) const -> VectorResult
{
    const unsigned R = m_reactions.size();
    const unsigned N = n.size();

    VectorResult res;
    func(res) = zeros(R);
    grad(res) = zeros(R, N);

    ScalarResult aux;
    for(unsigned i = 0; i < R; ++i)
    {
        aux = m_reactions[i].rate(T, P, n, a);

        func(res)[i]     = func(aux);
        grad(res).row(i) = grad(aux);
    }

    return res;
}

auto ReactionSystem::rates(const ChemicalState& state, const VectorResult& a) const -> VectorResult
{
    return rates(state.temperature(), state.pressure(), state.composition(), a);
}

auto ReactionSystem::initialiseStoichiometricMatrix(const ChemicalSystem& system) -> void
{
    auto species = system.species();

    const unsigned rows = m_reactions.size();
    const unsigned cols = system.numSpecies();

    stoichiometric_matrix = zeros(rows, cols);

    for(unsigned i = 0; i < rows; ++i)
        for(unsigned j = 0; j < cols; ++j)
            stoichiometric_matrix(i, j) = m_reactions[i].stoichiometry(species[j]);
}

auto ReactionSystem::initialiseReactionPhaseMap(const ChemicalSystem& system) -> void
{
    const Indices& species_phase_map = system.mapSpeciesToPhase();
    for(const Reaction& reaction : m_reactions)
    {
        std::set<Index> set_phase_indices;
        for(const Index& idx_species : reaction.idxReactingSpecies())
            set_phase_indices.insert(species_phase_map[idx_species]);
        Indices phase_indices(set_phase_indices.begin(), set_phase_indices.end());
        reaction_phase_map.push_back(phase_indices);
    }
}

auto ReactionSystem::initialiseReactionSpeciesMap() -> void
{
    for(const Reaction& reaction : m_reactions)
        reaction_species_map.push_back(reaction.idxReactingSpecies());
}

auto ReactionSystem::initialiseSpeciesReactionMap(const ChemicalSystem& system) -> void
{
    const unsigned num_species = system.numSpecies();
    const unsigned num_reactions = m_reactions.size();
    species_reaction_map.resize(num_species);
    for(unsigned i = 0; i < num_species; ++i)
        for(unsigned r = 0; r < num_reactions; ++r)
            if(contained(i, m_reactions[r].idxReactingSpecies()))
                species_reaction_map[i].push_back(r);
}

auto operator<<(std::ostream& out, const ReactionSystem& reactions) -> std::ostream&
{
    for(const Reaction& reaction : reactions.reactions())
        out << reaction << std::endl;

    return out;
}

} // namespace Reaktor
