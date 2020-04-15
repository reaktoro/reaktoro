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

#include "Reaction.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
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
    VectorXr stoichiometries;

    /// The function for the equilibrium constant of the reaction (in terms of natural log)
    std::function<real(real, real)> lnk;

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
        auto i = 0;
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

auto Reaction::setEquilibriumConstant(const std::function<real(real, real)>& lnk) -> void
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

auto Reaction::equilibriumConstant() const -> const std::function<real(real, real)>&
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

auto Reaction::stoichiometries() const -> VectorXrConstRef
{
    return pimpl->stoichiometries;
}

auto Reaction::stoichiometry(std::string species) const -> double
{
    return equation().stoichiometry(species);
}

auto Reaction::lnEquilibriumConstant(const ChemicalProps& props) const -> real
{
    // Get the temperature and pressure of the system
    const auto T = props.temperature();
    const auto P = props.pressure();

    // Check if a equilibrium constant function was provided
    if(pimpl->lnk) return pimpl->lnk(T, P);

    // Calculate the equilibrium constant using the standard Gibbs energies of the species
    const auto G0 = props.standardGibbsEnergies();
    const auto RT = universalGasConstant * T;

    real res = {};
    for(auto i = 0; i < indices().size(); ++i)
    {
        const auto ispecies = indices()[i];
        const auto vi = stoichiometries()[i];
        const auto G0i = G0[ispecies];
        res += vi * G0i;
    }

    return -res/RT;
}

auto Reaction::lnReactionQuotient(const ChemicalProps& props) const -> real
{
    const auto num_species = system().numSpecies();
    const auto ln_a = props.lnActivities();
    real ln_Q = {};
    auto counter = 0;
    for(auto i : indices())
    {
        const auto vi = stoichiometries()[counter];
        ln_Q += vi * ln_a[i];
        ++counter;
    }
    return ln_Q;
}

auto Reaction::lnEquilibriumIndex(const ChemicalProps& props) const -> real
{
	return lnReactionQuotient(props) - lnEquilibriumConstant(props);
}

auto Reaction::rate(const ChemicalProps& props) const -> real
{
    if(!pimpl->rate)
        errorFunctionNotInitialized("rate", "rate");
    return pimpl->rate(props);
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
