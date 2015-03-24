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
#include <Reaktor/Core/Element.hpp>
#include <Reaktor/Core/Utils.hpp>

namespace Reaktor {
namespace {

auto errorFunctionNotInitialized(std::string method) -> void
{
    std::string set_method = "set" + method + "Function";
    set_method[3] = std::toupper(set_method[3]);

    Exception exception;
    exception.error << "There was an error calling method `Species::" << method << "`.";
    exception.reason << "The error resulted because `Species::" << set_method << "` was not initialized called before.";
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
    ThermoScalarFunction standard_gibbs_energy_fn;

    /// The function for the apparent standard molar Helmholtz free energy of the species (in units of J/mol).
    ThermoScalarFunction standard_helmholtz_energy_fn;

    /// The function for the apparent standard molar internal energy of the species (in units of J/mol).
    ThermoScalarFunction standard_internal_energy_fn;

    /// The function for the apparent standard molar enthalpy of the species (in units of J/mol).
    ThermoScalarFunction standard_enthalpy_fn;

    /// The function for the standard molar entropy of the species (in units of J/K).
    ThermoScalarFunction standard_entropy_fn;

    /// The function for the standard molar volume of the species (in units of m3/mol).
    ThermoScalarFunction standard_volume_fn;

    /// The function for the standard molar isobaric heat capacity of the species (in units of J/(mol*K)).
    ThermoScalarFunction standard_heat_capacity_fn;
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

auto Species::setStandardGibbsEnergyFunction(const ThermoScalarFunction& function) -> void
{
    pimpl->standard_gibbs_energy_fn = function;
}

auto Species::setStandardHelmholtzEnergyFunction(const ThermoScalarFunction& function) -> void
{
    pimpl->standard_helmholtz_energy_fn = function;
}

auto Species::setStandardInternalEnergyFunction(const ThermoScalarFunction& function) -> void
{
    pimpl->standard_internal_energy_fn = function;
}

auto Species::setStandardEnthalpyFunction(const ThermoScalarFunction& function) -> void
{
    pimpl->standard_enthalpy_fn = function;
}

auto Species::setStandardEntropyFunction(const ThermoScalarFunction& function) -> void
{
    pimpl->standard_entropy_fn = function;
}

auto Species::setStandardVolumeFunction(const ThermoScalarFunction& function) -> void
{
    pimpl->standard_volume_fn = function;
}

auto Species::setStandardHeatCapacityFunction(const ThermoScalarFunction& function) -> void
{
    pimpl->standard_heat_capacity_fn = function;
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

auto Species::standardGibbsEnergyFunction() const -> const ThermoScalarFunction&
{
    return pimpl->standard_gibbs_energy_fn;
}

auto Species::standardHelmholtzEnergyFunction() const -> const ThermoScalarFunction&
{
    return pimpl->standard_helmholtz_energy_fn;
}

auto Species::standardInternalEnergyFunction() const -> const ThermoScalarFunction&
{
    return pimpl->standard_internal_energy_fn;
}

auto Species::standardEnthalpyFunction() const -> const ThermoScalarFunction&
{
    return pimpl->standard_enthalpy_fn;
}

auto Species::standardEntropyFunction() const -> const ThermoScalarFunction&
{
    return pimpl->standard_entropy_fn;
}

auto Species::standardVolumeFunction() const -> const ThermoScalarFunction&
{
    return pimpl->standard_volume_fn;
}

auto Species::standardHeatCapacityFunction() const -> const ThermoScalarFunction&
{
    return pimpl->standard_heat_capacity_fn;
}

auto Species::standardGibbsEnergy(double T, double P) const -> ThermoScalar
{
    if(not standardGibbsEnergyFunction())
        errorFunctionNotInitialized("standardGibbsEnergy");
    return standardGibbsEnergyFunction()(T, P);
}

auto Species::standardHelmholtzEnergy(double T, double P) const -> ThermoScalar
{
    if(not standardHelmholtzEnergyFunction())
        errorFunctionNotInitialized("standardHelmholtzEnergy");
    return standardHelmholtzEnergyFunction()(T, P);
}

auto Species::standardInternalEnergy(double T, double P) const -> ThermoScalar
{
    if(not standardInternalEnergyFunction())
        errorFunctionNotInitialized("standardInternalEnergy");
    return standardInternalEnergyFunction()(T, P);
}

auto Species::standardEnthalpy(double T, double P) const -> ThermoScalar
{
    if(not standardEnthalpyFunction())
        errorFunctionNotInitialized("standardEnthalpy");
    return standardEnthalpyFunction()(T, P);
}

auto Species::standardEntropy(double T, double P) const -> ThermoScalar
{
    if(not standardEntropyFunction())
        errorFunctionNotInitialized("standardEntropy");
    return standardEntropyFunction()(T, P);
}

auto Species::standardVolume(double T, double P) const -> ThermoScalar
{
    if(not standardVolumeFunction())
        errorFunctionNotInitialized("standardVolume");
    return standardVolumeFunction()(T, P);
}

auto Species::standardHeatCapacity(double T, double P) const -> ThermoScalar
{
    if(not standardHeatCapacityFunction())
        errorFunctionNotInitialized("standardHeatCapacity");
    return standardHeatCapacityFunction()(T, P);
}

auto operator<(const Species& lhs, const Species& rhs) -> bool
{
    return lhs.name() < rhs.name();
}

auto operator==(const Species& lhs, const Species& rhs) -> bool
{

    return lhs.name() == rhs.name();
}

} // namespace Reaktor
