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

#include "Reaction.hpp"

// Reaktoro includes
#include <Reaktoro/Common/ChemicalScalar.hpp>
#include <Reaktoro/Common/ChemicalVector.hpp>
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/ThermoScalar.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Core/ThermoProperties.hpp>

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

    /// The initial amounts of the species in the reaction (in units of mol)
    Vector n0;

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

auto Reaction::setInitialAmounts(VectorConstRef n_) -> void
{
    pimpl->n0 = n_;
}

auto Reaction::name() const -> std::string
{
    return pimpl->name;
}

auto Reaction::equilibriumConstant() const -> const ThermoScalarFunction&
{
    return pimpl->lnk;
}

auto Reaction::rate() const -> const ReactionRateFunction&
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

auto Reaction::initialAmounts() const -> const VectorRef
{
    return pimpl->n0;
}

auto Reaction::stoichiometries() const -> VectorConstRef
{
    return pimpl->stoichiometries;
}

auto Reaction::stoichiometry(std::string species) const -> double
{
    return equation().stoichiometry(species);
}

auto Reaction::lnEquilibriumConstant(const ChemicalProperties& properties) const -> ThermoScalar
{
    // Get the temperature and pressure of the system
    const double T = properties.temperature();
    const double P = properties.pressure();

    // Check if a equilibrium constant function was provided
    if(pimpl->lnk) return pimpl->lnk(T, P);

    // Calculate the equilibrium constant using the standard Gibbs energies of the species
    const ThermoVector G0 = properties.standardPartialMolarGibbsEnergies();
    const ThermoScalar RT = universalGasConstant * Temperature(T);

    ThermoScalar res;
    for(unsigned i = 0; i < indices().size(); ++i)
    {
        const Index ispecies = indices()[i];
        const double vi = stoichiometries()[i];
        const ThermoScalar G0i = G0[ispecies];
        res += vi * G0i;
    }

    return -res/RT;
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
        ln_Q += vi * ln_a[i];
        ++counter;
    }

    return ln_Q;
}

auto Reaction::lnEquilibriumIndex(const ChemicalProperties& properties) const -> ChemicalScalar
{
	return lnReactionQuotient(properties) - lnEquilibriumConstant(properties);
}

auto Reaction::rate(const ChemicalProperties& properties) const -> ChemicalScalar
{
    if(!pimpl->rate)
        errorFunctionNotInitialized("rate", "rate");
    return pimpl->rate(properties);
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
