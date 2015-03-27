// Reaktor is a C++ library for computational reaction modelling.
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

// Reaktor includes
#include <Reaktor/Common/Exception.hpp>
#include <Reaktor/Core/Multiphase.hpp>
#include <Reaktor/Core/Reaction.hpp>

namespace Reaktor {
namespace {

auto stoichiometricMatrix(const std::vector<Reaction>& reactions) -> Matrix
{
    const auto& multiphase = reactions.front().multiphase();
    const auto& species = multiphase.species();
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

    /// The multiphase instance
    Multiphase multiphase;

    /// The stoichiometric matrix of the reactions w.r.t. to all species in the system
    Matrix stoichiometric_matrix;

    /// The model configuration of the reaction system
    ReactionSystemModel model;

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

        // Initialize the multiphase instance
        multiphase = reactions.front().multiphase();

        // Initialize the stoichiometric matrix of the reactions
        stoichiometric_matrix = Reaktor::stoichiometricMatrix(reactions);

        const unsigned num_species = multiphase.numSpecies();
        const unsigned num_reactions = reactions.size();

        model.lnk = [&](double T, double P)
        {
            ThermoVector res(num_reactions);
            for(unsigned i = 0; i < num_reactions; ++i)
                res.row(i) = reactions[i].lnEquilibriumConstant(T, P);
            return res;
        };

        model.standard_gibbs_energy = [&](double T, double P)
        {
            ThermoVector res(num_reactions);
            for(unsigned i = 0; i < num_reactions; ++i)
                res.row(i) = reactions[i].standardGibbsEnergy(T, P);
            return res;
        };

        model.standard_helmholtz_energy = [&](double T, double P)
        {
            ThermoVector res(num_reactions);
            for(unsigned i = 0; i < num_reactions; ++i)
                res.row(i) = reactions[i].standardHelmholtzEnergy(T, P);
            return res;
        };

        model.standard_internal_energy = [&](double T, double P)
        {
            ThermoVector res(num_reactions);
            for(unsigned i = 0; i < num_reactions; ++i)
                res.row(i) = reactions[i].standardInternalEnergy(T, P);
            return res;
        };

        model.standard_enthalpy = [&](double T, double P)
        {
            ThermoVector res(num_reactions);
            for(unsigned i = 0; i < num_reactions; ++i)
                res.row(i) = reactions[i].standardEnthalpy(T, P);
            return res;
        };

        model.standard_entropy = [&](double T, double P)
        {
            ThermoVector res(num_reactions);
            for(unsigned i = 0; i < num_reactions; ++i)
                res.row(i) = reactions[i].standardEntropy(T, P);
            return res;
        };

        model.standard_volume = [&](double T, double P)
        {
            ThermoVector res(num_reactions);
            for(unsigned i = 0; i < num_reactions; ++i)
                res.row(i) = reactions[i].standardVolume(T, P);
            return res;
        };

        model.standard_heat_capacity = [&](double T, double P)
        {
            ThermoVector res(num_reactions);
            for(unsigned i = 0; i < num_reactions; ++i)
                res.row(i) = reactions[i].standardHeatCapacity(T, P);
            return res;
        };

        model.rates = [&](double T, double P, const Vector& n, const ChemicalVector& a)
        {
            ChemicalVector res(num_reactions, num_species);
            for(unsigned i = 0; i < num_reactions; ++i)
                res.row(i) = reactions[i].rate(T, P, n, a);
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

auto ReactionSystem::reactions() const -> const std::vector<Reaction>&
{
    return pimpl->reactions;
}

auto ReactionSystem::stoichiometricMatrix() const -> const Matrix&
{
    return pimpl->stoichiometric_matrix;
}

auto ReactionSystem::lnEquilibriumConstants(double T, double P) const -> ThermoVector
{
    return pimpl->model.lnk(T, P);
}

auto ReactionSystem::standardGibbsEnergies(double T, double P) const -> ThermoVector
{
    return pimpl->model.standard_gibbs_energy(T, P);
}

auto ReactionSystem::standardHelmholtzEnergies(double T, double P) const -> ThermoVector
{
    return pimpl->model.standard_helmholtz_energy(T, P);
}

auto ReactionSystem::standardInternalEnergies(double T, double P) const -> ThermoVector
{
    return pimpl->model.standard_internal_energy(T, P);
}

auto ReactionSystem::standardEnthalpies(double T, double P) const -> ThermoVector
{
    return pimpl->model.standard_enthalpy(T, P);
}

auto ReactionSystem::standardEntropies(double T, double P) const -> ThermoVector
{
    return pimpl->model.standard_entropy(T, P);
}

auto ReactionSystem::standardVolumes(double T, double P) const -> ThermoVector
{
    return pimpl->model.standard_volume(T, P);
}

auto ReactionSystem::standardHeatCapacities(double T, double P) const -> ThermoVector
{
    return pimpl->model.standard_heat_capacity(T, P);
}

auto ReactionSystem::rates(double T, double P, const Vector& n, const ChemicalVector& a) const -> ChemicalVector
{
    return pimpl->model.rates(T, P, n, a);
}

auto ReactionSystem::lnReactionQuotients(const ChemicalVector& a) const -> ChemicalVector
{
    const unsigned num_reactions = pimpl->reactions.size();
    const unsigned num_species = pimpl->multiphase.numSpecies();
    ChemicalVector res(num_reactions, num_species);
    for(unsigned i = 0; i < num_reactions; ++i)
        res.row(i) = reaction(i).lnReactionQuotient(a);
    return res;
}

} // namespace Reaktor
