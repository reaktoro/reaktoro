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

#include "Species.hpp"

// C++ includes
#include <set>

// Reaktor includes
#include <Reaktor/Common/Exception.hpp>
#include <Reaktor/Core/CoreUtils.hpp>
#include <Reaktor/Core/Element.hpp>

namespace Reaktor {
namespace {

auto errorFunctionNotInitialized(std::string method, std::string member) -> void
{
    Exception exception;
    exception.error << "There was an error calling method `Species::" << method << "`.";
    exception.reason << "The error resulted because `SpeciesData::" << member << "` was not initialized before constructing the Species instance.";
    RaiseError(exception);
}

} // namespace

struct Species::Impl
{
    /// The name of the chemical species
    std::string name;

    /// The chemical formula of the chemical species
    std::string formula;

    /// The elements that compose the chemical species and their coefficients
    std::map<Element, double> elements;

    /// The electrical charge of the chemical species
    double charge;

    /// The molar mass of the chemical species (in units of kg/mol)
    double molar_mass;

    /// The function for the apparent standard molar Gibbs free energy of the species (in units of J/mol).
    ThermoScalarFunction standard_gibbs_energy;

    /// The function for the apparent standard molar enthalpy of the species (in units of J/mol).
    ThermoScalarFunction standard_enthalpy;

    /// The function for the apparent standard molar Helmholtz free energy of the species (in units of J/mol).
    ThermoScalarFunction standard_helmholtz_energy;

    /// The function for the standard molar entropy of the species (in units of J/K).
    ThermoScalarFunction standard_entropy;

    /// The function for the standard molar volume of the species (in units of m3/mol).
    ThermoScalarFunction standard_volume;

    /// The function for the apparent standard molar internal energy of the species (in units of J/mol).
    ThermoScalarFunction standard_internal_energy;

    /// The function for the standard molar isobaric heat capacity of the species (in units of J/(mol*K)).
    ThermoScalarFunction standard_heat_capacity;
};

Species::Species()
: pimpl(new Impl())
{}

Species::Species(const Species& other)
: pimpl(new Impl(*other.pimpl))
{}

Species::~Species()
{}

auto Species::operator=(Species other) -> Species&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto Species::setName(std::string name) -> void
{
    pimpl->name = name;
}

auto Species::setFormula(std::string formula) -> void
{
    pimpl->formula = formula;
}

auto Species::setElements(const std::map<Element, double>& elements) -> void
{
    pimpl->elements = elements;
}

auto Species::setCharge(double value) -> void
{
    pimpl->charge = value;
}

auto Species::setMolarMass(double value) -> void
{
    pimpl->molar_mass = value;
}

auto Species::setStandardGibbsEnergy(const ThermoScalarFunction& function) -> void
{
    pimpl->standard_gibbs_energy = function;
}

auto Species::setStandardHelmholtzEnergy(const ThermoScalarFunction& function) -> void
{
    pimpl->standard_helmholtz_energy = function;
}

auto Species::setStandardInternalEnergy(const ThermoScalarFunction& function) -> void
{
    pimpl->standard_internal_energy = function;
}

auto Species::setStandardEnthalpy(const ThermoScalarFunction& function) -> void
{
    pimpl->standard_enthalpy = function;
}

auto Species::setStandardEntropy(const ThermoScalarFunction& function) -> void
{
    pimpl->standard_entropy = function;
}

auto Species::setStandardVolume(const ThermoScalarFunction& function) -> void
{
    pimpl->standard_volume = function;
}

auto Species::setStandardHeatCapacity(const ThermoScalarFunction& function) -> void
{
    pimpl->standard_heat_capacity = function;
}

auto Species::numElements() const -> unsigned
{
    return elements().size();
}

auto Species::name() const -> const std::string&
{
    return pimpl->name;
}

auto Species::formula() const -> const std::string&
{
    return pimpl->formula;
}

auto Species::elements() const -> const std::map<Element, double>&
{
    return pimpl->elements;
}

auto Species::charge() const -> double
{
    return pimpl->charge;
}

auto Species::molarMass() const -> double
{
    return pimpl->molar_mass;
}

auto Species::elementAtoms(std::string element) const -> double
{
    for(const auto& pair : elements())
        if(element == pair.first.name())
            return pair.second;
    return 0.0;
}

auto Species::standardGibbsEnergy(double T, double P) const -> ThermoScalar
{
    if(not pimpl->standard_gibbs_energy)
        errorFunctionNotInitialized("standardGibbsEnergy", "standard_gibbs_energy");
    return pimpl->standard_gibbs_energy(T, P);
}

auto Species::standardHelmholtzEnergy(double T, double P) const -> ThermoScalar
{
    if(not pimpl->standard_helmholtz_energy)
        errorFunctionNotInitialized("standardHelmholtzEnergy", "standard_helmholtz_energy");
    return pimpl->standard_helmholtz_energy(T, P);
}

auto Species::standardInternalEnergy(double T, double P) const -> ThermoScalar
{
    if(not pimpl->standard_internal_energy)
        errorFunctionNotInitialized("standardInternalEnergy", "standard_internal_energy");
    return pimpl->standard_internal_energy(T, P);
}

auto Species::standardEnthalpy(double T, double P) const -> ThermoScalar
{
    if(not pimpl->standard_enthalpy)
        errorFunctionNotInitialized("standardEnthalpy", "standard_enthalpy");
    return pimpl->standard_enthalpy(T, P);
}

auto Species::standardEntropy(double T, double P) const -> ThermoScalar
{
    if(not pimpl->standard_entropy)
        errorFunctionNotInitialized("standardEntropy", "standard_entropy");
    return pimpl->standard_entropy(T, P);
}

auto Species::standardVolume(double T, double P) const -> ThermoScalar
{
    if(not pimpl->standard_volume)
        errorFunctionNotInitialized("standardVolume", "standard_volume");
    return pimpl->standard_volume(T, P);
}

auto Species::standardHeatCapacity(double T, double P) const -> ThermoScalar
{
    if(not pimpl->standard_heat_capacity)
        errorFunctionNotInitialized("standardHeatCapacity", "standard_heat_capacity");
    return pimpl->standard_heat_capacity(T, P);
}

auto operator<(const Species& lhs, const Species& rhs) -> bool
{
    return lhs.name() < rhs.name();
}

auto operator==(const Species& lhs, const Species& rhs) -> bool
{

    return lhs.name() == rhs.name();
}

auto collectElements(const std::vector<Species>& species) -> std::vector<Element>
{
    std::set<Element> elements;
    for(const Species& iter : species)
        for(const auto& pair : iter.elements())
            elements.insert(pair.first);
    return std::vector<Element>(elements.begin(), elements.end());
}

auto charges(const std::vector<Species>& species) -> Vector
{
    Vector charges(species.size());
    for(unsigned i = 0; i < species.size(); ++i)
        charges[i] = species[i].charge();
    return charges;
}

auto molarMasses(const std::vector<Species>& species) -> Vector
{
    Vector molar_masses(species.size());
    for(unsigned i = 0; i < species.size(); ++i)
        molar_masses[i] = species[i].molarMass();
    return molar_masses;
}

} // namespace Reaktor
