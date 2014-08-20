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

#include "Reaction.hpp"

// C++ includes
#include <cmath>
#include <iomanip>
#include <sstream>

// Reaktor includes
#include <Reaktor/Common/Macros.hpp>
#include <Reaktor/Core/ChemicalSystem.hpp>
#include <Reaktor/Core/ReactionEquation.hpp>
#include <Reaktor/Core/Species.hpp>
#include <Reaktor/Utils/SetUtils.hpp>

namespace Reaktor {
namespace internal {

/**
 * Creates a default equilibrium constant function for a reaction
 *
 * The created function will use the chemical potential functions of the
 * participating species to calculate the equilibrium constant of the reaction.
 */
auto defaultEquilibriumConstant(const ReactionEquation& equation, const ChemicalSystem& system) -> EquilibriumConstantFn
{
    // The chemical species in the system
    const std::vector<Species>& species = system.species();

    // The chemical potentials of the species in the reaction
    std::vector<ChemicalPotentialFn> chemical_potentials;

    // The stoichiometries of the species in the reaction
    std::vector<double> stoichiometries;

    for(const auto& pair : equation)
    {
        const Index ispecies = system.idxSpeciesWithError(pair.first);
        chemical_potentials.push_back(species[ispecies].chemicalPotential());
        stoichiometries.push_back(pair.second);
    }

    // The number of participating species in the reaction
    const unsigned num_species = equation.size();

    // The universal gas constant (in units of J/(mol*K))
    const double R = 8.3144621;

    // Define the equilibrium constant function
    EquilibriumConstantFn K = [=](double T, double P)
    {
        double sum = 0.0;
        for(unsigned i = 0; i < num_species; ++i)
            sum += stoichiometries[i] * chemical_potentials[i](T, P);

        return std::exp(-sum/(R*T));
    };

    return K;
}

/**
 * Creates a default (zero) reaction rate function for a reaction
 */
auto defaultRate(const ChemicalSystem& system) -> ReactionRateFn
{
    const Vector zero = zeros(system.numSpecies());
    const PartialScalar rate = partialScalar(0.0, zero);

    ReactionRateFn fn = [=](double T, double P, const Vector& n, const PartialVector& a)
    {
        return rate;
    };

    return fn;
}

} /* namespace internal */

using namespace internal;

class Reaction::Impl
{
public:
    /// The reaction equation
    ReactionEquation equation;

    /// The species that participate in the reaction
    std::vector<std::string> species;

    /// The stoichiometries of the species
    std::vector<double> stoichiometries;

    /// The indices of the species in the reaction
    Indices idx_species;

    /// The equilibrium constant of the reaction
    EquilibriumConstantFn equilibrium_constant;

    /// The rate function of the reaction
    ReactionRateFn rate;

public:
    Impl()
    {}

    Impl(const ReactionEquation& equation, const ChemicalSystem& system)
    : equation(equation)
    {
        // Initialise the species, their indices and stoichiometries
        for(const auto& pair : equation)
        {
            species.push_back(pair.first);
            idx_species.push_back(system.idxSpeciesWithError(pair.first));
            stoichiometries.push_back(pair.second);
        }

        // Set the default equilibrium constant of the reaction
        setEquilibriumConstantFn(defaultEquilibriumConstant(equation, system));

        // Set the default rate of the reaction
        setRate(defaultRate(system));
    }

    auto setEquilibriumConstantFn(const EquilibriumConstantFn& equilibrium_constant) -> void
    {
        this->equilibrium_constant = equilibrium_constant;
    }

    auto setRate(const ReactionRateFn& rate) -> void
    {
        this->rate = rate;
    }

    auto stoichiometry(const std::string& name) const -> double
    {
        const unsigned idx = find(name, species);
        return (idx < species.size()) ? stoichiometries[idx] : 0.0;
    }

    auto equilibriumConstant(double T, double P) const -> double
    {
        return equilibrium_constant(T, P);
    }

    auto reactionQuotient(const PartialVector& a) const -> PartialScalar
    {
        const unsigned N = func(a).size();

        PartialScalar Q = partialScalar(1.0, zeros(N));

        for(unsigned i = 0; i < idx_species.size(); ++i)
        {
            const double vi = stoichiometries[i];
            const double ai = func(a)[idx_species[i]];
            func(Q) *= std::pow(ai, vi);
        }

        for(unsigned i = 0; i < idx_species.size(); ++i)
        {
            const double vi = stoichiometries[i];
            const double ai = func(a)[idx_species[i]];
            grad(Q) += func(Q) * vi/ai * grad(a).row(idx_species[i]);
        }

        return Q;
    }
};

Reaction::Reaction()
: pimpl(new Impl())
{}

Reaction::Reaction(const ReactionEquation& equation, const ChemicalSystem& system)
: pimpl(new Impl(equation, system))
{}

Reaction::Reaction(const Reaction& other)
: pimpl(new Impl(*other.pimpl))
{}

Reaction::~Reaction()
{}

auto Reaction::operator=(Reaction other) -> Reaction&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto Reaction::setEquilibriumConstant(const EquilibriumConstantFn& equilibrium_constant) -> void
{
    pimpl->setEquilibriumConstantFn(equilibrium_constant);
}

auto Reaction::setRate(const ReactionRateFn& rate) -> void
{
    pimpl->setRate(rate);
}

auto Reaction::equation() const -> const ReactionEquation&
{
    return pimpl->equation;
}

auto Reaction::rate() const -> const ReactionRateFn&
{
    return pimpl->rate;
}

auto Reaction::equilibriumConstant() const -> const EquilibriumConstantFn&
{
    return pimpl->equilibrium_constant;
}

auto Reaction::idxReactingSpecies() const -> const Indices&
{
    return pimpl->idx_species;
}

auto Reaction::stoichiometry(const std::string& species) const -> double
{
    return pimpl->stoichiometry(species);
}

auto Reaction::equilibriumConstant(double T, double P) const -> double
{
    return pimpl->equilibriumConstant(T, P);
}

auto Reaction::reactionQuotient(const PartialVector& a) const -> PartialScalar
{
    return pimpl->reactionQuotient(a);
}

auto Reaction::rate(double T, double P, const Vector& n, const PartialVector& a) const -> PartialScalar
{
    return pimpl->rate(T, P, n, a);
}

auto operator<<(std::ostream& out, const Reaction& reaction) -> std::ostream&
{
    out << reaction.equation();

    return out;
}

} /* namespace Reaktor */
