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

#pragma once

// C++ includes
#include <map>
#include <memory>
#include <string>
#include <vector>

// Reaktor includes
#include <Reaktor/Common/Matrix.hpp>
#include <Reaktor/Common/ThermoScalar.hpp>

namespace Reaktor {

// Forward declarations
class Element;

/// A type used to describe a chemical species and its attributes.
/// The Species class is used to represent a chemical species. It is an important
/// class in the library, since it defines fundamental attributes of a general
/// chemical species such as its elemental formula, electrical charge and molar mass.
/// @see Phase
/// @ingroup Core
class Species
{
public:
    /// Construct a default Species instance.
    Species();

    /// Construct a copy of an Species instance
    Species(const Species& other);

    /// Destroy this instance
    virtual ~Species();

    /// Assign an Species instance to this instance
    auto operator=(Species other) -> Species&;

    /// Set the name of the species.
    auto setName(std::string name) -> void;

    /// Set the formula of the species.
    auto setFormula(std::string formula) -> void;

    /// Set the elements of the species.
    auto setElements(const std::map<Element, double>& elements) -> void;

    /// Set the charge of the species.
    auto setCharge(double value) -> void;

    /// Set the molar mass of the species (in units of kg/mol).
    auto setMolarMass(double value) -> void;

    /// Set the standard Gibbs energy function of the species.
    auto setStandardGibbsEnergyFunction(const ThermoScalarFunction& function) -> void;

    /// Set the standard Helmholtz energy function of the species.
    auto setStandardHelmholtzEnergyFunction(const ThermoScalarFunction& function) -> void;

    /// Set the standard internal energy function of the species.
    auto setStandardInternalEnergyFunction(const ThermoScalarFunction& function) -> void;

    /// Set the standard enthalpy function of the species.
    auto setStandardEnthalpyFunction(const ThermoScalarFunction& function) -> void;

    /// Set the standard entropy function of the species.
    auto setStandardEntropyFunction(const ThermoScalarFunction& function) -> void;

    /// Set the standard volume function of the species.
    auto setStandardVolumeFunction(const ThermoScalarFunction& function) -> void;

    /// Set the standard heat capacity function of the species.
    auto setStandardHeatCapacityFunction(const ThermoScalarFunction& function) -> void;

    /// Return the number of elements of the chemical species
    auto numElements() const -> unsigned;

    /// Return the name of the chemical species
    auto name() const -> const std::string&;

    /// Return the formula of the chemical species
    auto formula() const -> const std::string&;

    /// Return the elements that compose the chemical species and their coefficients
    auto elements() const -> const std::map<Element, double>&;

    /// Return the electrical charge of the chemical species
    auto charge() const -> double;

    /// Return the molar mass of the chemical species (in units of kg/mol)
    auto molarMass() const -> double;

    /// Return the number of atoms of an element in the chemical species.
    auto elementAtoms(std::string element) const -> double;

    /// Return the standard Gibbs energy function of the species.
    auto standardGibbsEnergyFunction() const -> const ThermoScalarFunction&;

    /// Return the standard Helmholtz energy function of the species.
    auto standardHelmholtzEnergyFunction() const -> const ThermoScalarFunction&;

    /// Return the standard internal energy function of the species.
    auto standardInternalEnergyFunction() const -> const ThermoScalarFunction&;

    /// Return the standard enthalpy function of the species.
    auto standardEnthalpyFunction() const -> const ThermoScalarFunction&;

    /// Return the standard entropy function of the species.
    auto standardEntropyFunction() const -> const ThermoScalarFunction&;

    /// Return the standard volume function of the species.
    auto standardVolumeFunction() const -> const ThermoScalarFunction&;

    /// Return the standard heat capacity function of the species.
    auto standardHeatCapacityFunction() const -> const ThermoScalarFunction&;

    /// Calculate the apparent standard molar Gibbs free energy of the species (in units of J/mol).
    auto standardGibbsEnergy(double T, double P) const -> ThermoScalar;

    /// Calculate the apparent standard molar Helmholtz free energy of the species (in units of J/mol).
    auto standardHelmholtzEnergy(double T, double P) const -> ThermoScalar;

    /// Calculate the apparent standard molar internal energy of the species (in units of J/mol).
    auto standardInternalEnergy(double T, double P) const -> ThermoScalar;

    /// Calculate the apparent standard molar enthalpy of the species (in units of J/mol).
    auto standardEnthalpy(double T, double P) const -> ThermoScalar;

    /// Calculate the standard molar entropies of the species (in units of J/K).
    auto standardEntropy(double T, double P) const -> ThermoScalar;

    /// Calculate the standard molar volumes of the species (in units of m3/mol).
    auto standardVolume(double T, double P) const -> ThermoScalar;

    /// Calculate the standard molar isobaric heat capacity of the species (in units of J/(mol*K)).
    auto standardHeatCapacity(double T, double P) const -> ThermoScalar;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

/// Compare two Species instances for less than
auto operator<(const Species& lhs, const Species& rhs) -> bool;

/// Compare two Species instances for equality
auto operator==(const Species& lhs, const Species& rhs) -> bool;

} // namespace Reaktor
