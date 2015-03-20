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

#include "Reaction.hpp"

// Reaktor includes
#include <Reaktor/Common/ChemicalScalar.hpp>
#include <Reaktor/Common/ChemicalVector.hpp>
#include <Reaktor/Core/ChemicalSystem.hpp>
#include <Reaktor/Common/Exception.hpp>

namespace Reaktor {
namespace {

auto errorFunctionNotInitialized(std::string method, std::string member) -> void
{
    Exception exception;
    exception.error << "There was an error calling method `Reaktor::" << method << "`.";
    exception.reason << "The error resulted because `ReactionData::" << member << "` was not initialized before constructing the Reaction instance.";
    RaiseError(exception);
}

} // namespace

struct Reaction::Impl
{
    /// The data of the reaction
    ReactionData data;

    /// The names of the reacting species of the reaction
    std::vector<std::string> species;

    /// The indices of the reacting species of the reaction
    Indices indices;

    /// The stoichiometries of the reacting species of the reaction
    std::vector<double> stoichiometries;
};

Reaction::Reaction()
: pimpl(new Impl())
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

auto Reaction::equation() const -> const ReactionEquation&
{
    return pimpl->data.equation;
}

auto Reaction::species() const -> const std::vector<std::string>&
{
    return pimpl->species;
}

auto Reaction::indices() const -> const Indices&
{
    return pimpl->indices;
}

auto Reaction::stoichiometries() const -> const std::vector<double>&
{
    return pimpl->stoichiometries;
}

auto Reaction::lnEquilibriumConstant(double T, double P) const -> ThermoScalar
{
    if(not pimpl->data.lnk)
        errorFunctionNotInitialized("lnEquilibriumConstant", "lnk");
    return pimpl->data.lnk(T, P);
}

auto Reaction::standardGibbsEnergy(double T, double P) const -> ThermoScalar
{
    if(not pimpl->data.standard_gibbs_energy)
        errorFunctionNotInitialized("standardGibbsEnergy", "standard_gibbs_energy");
    return pimpl->data.standard_gibbs_energy(T, P);
}

auto Reaction::standardHelmholtzEnergy(double T, double P) const -> ThermoScalar
{
    if(not pimpl->data.standard_helmholtz_energy)
        errorFunctionNotInitialized("standardHelmholtzEnergy", "standard_helmholtz_energy");
    return pimpl->data.standard_helmholtz_energy(T, P);
}

auto Reaction::standardInternalEnergy(double T, double P) const -> ThermoScalar
{
    if(not pimpl->data.standard_internal_energy)
        errorFunctionNotInitialized("standardInternalEnergy", "standard_internal_energy");
    return pimpl->data.standard_internal_energy(T, P);
}

auto Reaction::standardEnthalpy(double T, double P) const -> ThermoScalar
{
    if(not pimpl->data.standard_enthalpy)
        errorFunctionNotInitialized("standardEnthalpy", "standard_enthalpy");
    return pimpl->data.standard_enthalpy(T, P);
}

auto Reaction::standardEntropy(double T, double P) const -> ThermoScalar
{
    if(not pimpl->data.standard_entropy)
        errorFunctionNotInitialized("standardEntropy", "standard_entropy");
    return pimpl->data.standard_entropy(T, P);
}

auto Reaction::standardVolume(double T, double P) const -> ThermoScalar
{
    if(not pimpl->data.standard_volume)
        errorFunctionNotInitialized("standardVolume", "standard_volume");
    return pimpl->data.standard_volume(T, P);
}

auto Reaction::standardHeatCapacity(double T, double P) const -> ThermoScalar
{
    if(not pimpl->data.standard_heat_capacity)
        errorFunctionNotInitialized("standardHeatCapacity", "standard_heat_capacity");
    return pimpl->data.standard_heat_capacity(T, P);
}

auto Reaction::rate(double T, double P, const Vector& n, const ChemicalVector& ln_a) const -> ChemicalScalar
{
    if(not pimpl->data.rate)
        errorFunctionNotInitialized("rate", "rate");
    return pimpl->data.rate(T, P, n, ln_a);
}

auto Reaction::lnReactionQuotient(const ChemicalVector& ln_a) const -> ChemicalScalar
{
    const unsigned num_species = ln_a.val.size();
    ChemicalScalar lnQ(num_species);

    unsigned counter = 0;
    for(Index i : indices())
    {
        const double vi = stoichiometries()[counter];
        lnQ.val += vi * ln_a.val[i];
        lnQ.ddn += vi * ln_a.ddn.row(i);
        ++counter;
    }

    return lnQ;
}

} // namespace Reaktor
