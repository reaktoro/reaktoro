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
#include <Reaktor/Common/ChemicalVector.hpp>
#include <Reaktor/Common/Constants.hpp>
#include <Reaktor/Common/Exception.hpp>
#include <Reaktor/Core/ChemicalSystem.hpp>

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
    /// The equation of the reaction as a list of species and stoichiometries
    ReactionEquation equation;

    /// The chemical system instance
    ChemicalSystem system;

    /// The species in the reaction
    std::vector<Species> species;

    /// The indices of the species in the reaction
    Indices indices;

    /// The stoichiometries of the species in the reaction
    Vector stoichiometries;

    /// The function for the equilibrium constant of the reaction (in terms of natural log)
    ThermoScalarFunction lnk;

    /// The function for the apparent standard molar Gibbs free energy of the reaction (in units of J/mol).
    ThermoScalarFunction standard_gibbs_energy;

    /// The function for the apparent standard molar Helmholtz free energy of the reaction (in units of J/mol).
    ThermoScalarFunction standard_helmholtz_energy;

    /// The function for the apparent standard molar internal energy of the reaction (in units of J/mol).
    ThermoScalarFunction standard_internal_energy;

    /// The function for the apparent standard molar enthalpy of the reaction (in units of J/mol).
    ThermoScalarFunction standard_enthalpy;

    /// The function for the standard molar entropy of the reaction (in units of J/K).
    ThermoScalarFunction standard_entropy;

    /// The function for the standard molar volume of the reaction (in units of m3/mol).
    ThermoScalarFunction standard_volume;

    /// The function for the standard molar isobaric heat capacity of the reaction (in units of J/(mol*K)).
    ThermoScalarFunction standard_heat_capacity;

    /// The function for the kinetic rate of the reaction (in units of mol/s)
    ReactionRateFunction rate;

    Impl()
    {}

    Impl(const ReactionEquation& equation, const ChemicalSystem& system)
    : equation(equation), system(system)
    {
        // Initialize the species, their indices, and their stoichiometries in the reaction
        species.resize(equation.size());
        indices.resize(equation.size());
        stoichiometries.resize(equation.size());
        for(unsigned i = 0; i < equation.size(); ++i)
        {
            species[i] = system.species(equation[i].first);
            indices[i] = system.indexSpecies(equation[i].first);
            stoichiometries[i] = equation[i].second;
        }

        // Initialize the function for the apparent standard molar Gibbs free energy of the reaction
        standard_gibbs_energy = [=](double T, double P) -> ThermoScalar
        {
            ThermoScalar res;
            for(unsigned i = 0; i < indices.size(); ++i)
                res += stoichiometries[i] * species[i].standardGibbsEnergy(T, P);
            return res;
        };

        // Initialize the function for the apparent standard molar Helmholtz free energy of the reaction
        standard_helmholtz_energy = [=](double T, double P) -> ThermoScalar
        {
            ThermoScalar res;
            for(unsigned i = 0; i < indices.size(); ++i)
                res += stoichiometries[i] * species[i].standardHelmholtzEnergy(T, P);
            return res;
        };

        // Initialize the function for the apparent standard molar internal energy of the reaction
        standard_internal_energy = [=](double T, double P) -> ThermoScalar
        {
            ThermoScalar res;
            for(unsigned i = 0; i < indices.size(); ++i)
                res += stoichiometries[i] * species[i].standardInternalEnergy(T, P);
            return res;
        };

        // Initialize the function for the apparent standard molar enthalpy of the reaction
        standard_enthalpy = [=](double T, double P) -> ThermoScalar
        {
            ThermoScalar res;
            for(unsigned i = 0; i < indices.size(); ++i)
                res += stoichiometries[i] * species[i].standardEnthalpy(T, P);
            return res;
        };

        // Initialize the function for the standard molar entropy of the reaction
        standard_entropy = [=](double T, double P) -> ThermoScalar
        {
            ThermoScalar res;
            for(unsigned i = 0; i < indices.size(); ++i)
                res += stoichiometries[i] * species[i].standardEntropy(T, P);
            return res;
        };

        // Initialize the function for the standard molar volume of the reaction
        standard_volume = [=](double T, double P) -> ThermoScalar
        {
            ThermoScalar res;
            for(unsigned i = 0; i < indices.size(); ++i)
                res += stoichiometries[i] * species[i].standardVolume(T, P);
            return res;
        };

        // Initialize the function for the standard molar isobaric heat capacity of the reaction
        standard_heat_capacity = [=](double T, double P) -> ThermoScalar
        {
            ThermoScalar res;
            for(unsigned i = 0; i < indices.size(); ++i)
                res += stoichiometries[i] * species[i].standardHeatCapacity(T, P);
            return res;
        };

        // Initialize the function for the equilibrium constant of the reaction
        lnk = [=](double T, double P) -> ThermoScalar
        {
            const double R = universalGasConstant;
            return -standard_gibbs_energy(T, P)/(R*T);
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

auto Reaction::setEquilibriumConstantFunction(const ThermoScalarFunction& lnk) -> void
{
    pimpl->lnk = lnk;
}

auto Reaction::setStandardGibbsEnergyFunction(const ThermoScalarFunction& function) -> void
{
    pimpl->standard_gibbs_energy = function;
}

auto Reaction::setStandardHelmholtzEnergyFunction(const ThermoScalarFunction& function) -> void
{
    pimpl->standard_helmholtz_energy = function;
}

auto Reaction::setStandardInternalEnergyFunction(const ThermoScalarFunction& function) -> void
{
    pimpl->standard_internal_energy = function;
}

auto Reaction::setStandardEnthalpyFunction(const ThermoScalarFunction& function) -> void
{
    pimpl->standard_enthalpy = function;
}

auto Reaction::setStandardEntropyFunction(const ThermoScalarFunction& function) -> void
{
    pimpl->standard_entropy = function;
}

auto Reaction::setStandardVolumeFunction(const ThermoScalarFunction& function) -> void
{
    pimpl->standard_volume = function;
}

auto Reaction::setStandardHeatCapacityFunction(const ThermoScalarFunction& function) -> void
{
    pimpl->standard_heat_capacity = function;
}

auto Reaction::setRate(const ReactionRateFunction& function) -> void
{
    pimpl->rate = function;
}

auto Reaction::equilibriumConstantFunction() const -> ThermoScalarFunction
{
    return pimpl->lnk;
}

auto Reaction::standardGibbsEnergyFunction() const -> ThermoScalarFunction
{
    return pimpl->standard_gibbs_energy;
}

auto Reaction::standardHelmholtzEnergyFunction() const -> ThermoScalarFunction
{
    return pimpl->standard_helmholtz_energy;
}

auto Reaction::standardInternalEnergyFunction() const -> ThermoScalarFunction
{
    return pimpl->standard_internal_energy;
}

auto Reaction::standardEnthalpyFunction() const -> ThermoScalarFunction
{
    return pimpl->standard_enthalpy;
}

auto Reaction::standardEntropyFunction() const -> ThermoScalarFunction
{
    return pimpl->standard_entropy;
}

auto Reaction::standardVolumeFunction() const -> ThermoScalarFunction
{
    return pimpl->standard_volume;
}

auto Reaction::standardHeatCapacityFunction() const -> ThermoScalarFunction
{
    return pimpl->standard_heat_capacity;
}

auto Reaction::rateFunction() const -> ReactionRateFunction
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
    if(not pimpl->lnk)
        errorFunctionNotInitialized("lnEquilibriumConstant", "lnk");
    return pimpl->lnk(T, P);
}

auto Reaction::standardGibbsEnergy(double T, double P) const -> ThermoScalar
{
    if(not pimpl->standard_gibbs_energy)
        errorFunctionNotInitialized("standardGibbsEnergy", "standard_gibbs_energy");
    return pimpl->standard_gibbs_energy(T, P);
}

auto Reaction::standardHelmholtzEnergy(double T, double P) const -> ThermoScalar
{
    if(not pimpl->standard_helmholtz_energy)
        errorFunctionNotInitialized("standardHelmholtzEnergy", "standard_helmholtz_energy");
    return pimpl->standard_helmholtz_energy(T, P);
}

auto Reaction::standardInternalEnergy(double T, double P) const -> ThermoScalar
{
    if(not pimpl->standard_internal_energy)
        errorFunctionNotInitialized("standardInternalEnergy", "standard_internal_energy");
    return pimpl->standard_internal_energy(T, P);
}

auto Reaction::standardEnthalpy(double T, double P) const -> ThermoScalar
{
    if(not pimpl->standard_enthalpy)
        errorFunctionNotInitialized("standardEnthalpy", "standard_enthalpy");
    return pimpl->standard_enthalpy(T, P);
}

auto Reaction::standardEntropy(double T, double P) const -> ThermoScalar
{
    if(not pimpl->standard_entropy)
        errorFunctionNotInitialized("standardEntropy", "standard_entropy");
    return pimpl->standard_entropy(T, P);
}

auto Reaction::standardVolume(double T, double P) const -> ThermoScalar
{
    if(not pimpl->standard_volume)
        errorFunctionNotInitialized("standardVolume", "standard_volume");
    return pimpl->standard_volume(T, P);
}

auto Reaction::standardHeatCapacity(double T, double P) const -> ThermoScalar
{
    if(not pimpl->standard_heat_capacity)
        errorFunctionNotInitialized("standardHeatCapacity", "standard_heat_capacity");
    return pimpl->standard_heat_capacity(T, P);
}

auto Reaction::rate(double T, double P, const Vector& n, const ChemicalVector& a) const -> ChemicalScalar
{
    if(not pimpl->rate)
        errorFunctionNotInitialized("rate", "rate");
    return pimpl->rate(T, P, n, a);
}

auto Reaction::lnReactionQuotient(const ChemicalVector& a) const -> ChemicalScalar
{
    const unsigned num_species = system().numSpecies();
    ChemicalVector lna = log(a);
    ChemicalScalar lnQ(num_species);

    unsigned counter = 0;
    for(Index i : indices())
    {
        const double vi = stoichiometries()[counter];
        lnQ += vi * lna.row(i);
        ++counter;
    }

    return lnQ;
}

} // namespace Reaktor
