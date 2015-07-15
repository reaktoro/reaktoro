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

#include "Reaction.hpp"

// Reaktoro includes
#include <Reaktoro/Common/ChemicalVector.hpp>
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Core/SpeciesProperties.hpp>

namespace Reaktoro {
namespace {

auto errorFunctionNotInitialized(std::string method, std::string member) -> void
{
    Exception exception;
    exception.error << "There was an error calling method `Reaktoro::" << method << "`.";
    exception.reason << "The error resulted because `ReactionData::" << member << "` was not initialized before constructing the Reaction instance.";
    RaiseError(exception);
}

} // namespace

struct Reaction::Impl
{
    /// The equation of the reaction as a list of species and stoichiometries
    ReactionEquation equation;

    /// The chemical system instance
    ChemicalSystem system;

    /// The name of the reaction
    std::string name;

    /// The species in the reaction
    std::vector<Species> species;

    /// The indices of the species in the reaction
    Indices indices;

    /// The stoichiometries of the species in the reaction
    Vector stoichiometries;

    /// The function for the equilibrium constant of the reaction (in terms of natural log)
    ThermoScalarFunction lnk;

    /// The function for the kinetic rate of the reaction (in units of mol/s)
    ReactionRateFunction rate;

    Impl()
    {}

    Impl(const ReactionEquation& equation, const ChemicalSystem& system)
    : equation(equation), system(system)
    {
        // Initialize the species, their indices, and their stoichiometries in the reaction
        species.reserve(equation.numSpecies());
        indices.reserve(equation.numSpecies());
        stoichiometries.resize(equation.numSpecies());
        unsigned i = 0;
        for(const auto& pair : equation.equation())
        {
            species.push_back(system.species(pair.first));
            indices.push_back(system.indexSpecies(pair.first));
            stoichiometries[i] = pair.second;
            ++i;
        }

        // Define auxiliary copies to avoid reference to this pointer in the lambdas below
        const auto species_ = species;
        const auto stoichiometries_ = stoichiometries;

        // Initialize the function for the equilibrium constant of the reaction
        lnk = [=](double T, double P) -> ThermoScalar
        {
            ThermoScalar res;
            for(unsigned i = 0; i < species_.size(); ++i)
            {
                SpeciesProperties sp = species_[i].properties(T, P);
                res += stoichiometries_[i] * sp.standardPartialMolarGibbsEnergy();
            }

            const ThermoScalar RT = universalGasConstant*ThermoScalar(T, 1.0, 0.0);
            return -res/RT;
        };
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

auto Reaction::setName(std::string name) -> void
{
    pimpl->name = name;
}

auto Reaction::setEquilibriumConstant(const ThermoScalarFunction& lnk) -> void
{
    pimpl->lnk = lnk;
}

auto Reaction::setRate(const ReactionRateFunction& function) -> void
{
    pimpl->rate = function;
}

auto Reaction::name() const -> std::string
{
    return pimpl->name;
}

auto Reaction::equilibriumConstantFunction() const -> const ThermoScalarFunction&
{
    return pimpl->lnk;
}

auto Reaction::rateFunction() const -> const ReactionRateFunction&
{
    return pimpl->rate;
}

auto Reaction::equation() const -> const ReactionEquation&
{
    return pimpl->equation;
}

auto Reaction::system() const -> const ChemicalSystem&
{
    return pimpl->system;
}

auto Reaction::species() const -> const std::vector<Species>&
{
    return pimpl->species;
}

auto Reaction::indices() const -> const Indices&
{
    return pimpl->indices;
}

auto Reaction::stoichiometries() const -> const Vector&
{
    return pimpl->stoichiometries;
}

auto Reaction::stoichiometry(std::string species) const -> double
{
    return equation().stoichiometry(species);
}

auto Reaction::lnEquilibriumConstant(double T, double P) const -> ThermoScalar
{
    if(!pimpl->lnk)
        errorFunctionNotInitialized("lnEquilibriumConstant", "lnk");
    return pimpl->lnk(T, P);
}

auto Reaction::rate(const ChemicalProperties& properties) const -> ChemicalScalar
{
    if(!pimpl->rate)
        errorFunctionNotInitialized("rate", "rate");
    return pimpl->rate(properties);
}

auto Reaction::lnReactionQuotient(const ChemicalProperties& properties) const -> ChemicalScalar
{
    const unsigned num_species = system().numSpecies();
    const ChemicalVector& ln_a = properties.lnActivities();
    ChemicalScalar ln_Q(num_species);
    unsigned counter = 0;
    for(Index i : indices())
    {
        const double vi = stoichiometries()[counter];
        ln_Q += vi * ln_a.row(i);
        ++counter;
    }

    return ln_Q;
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
