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
#include <memory>
#include <string>
#include <vector>

// Reaktor includes
#include <Reaktor/Common/ThermoProperty.hpp>

namespace Reaktor {

// Forward declarations
struct SpeciesThermoModel;

/// A type used to describe a chemical species and its attributes.
/// The Species class is used to represent a chemical species. It is an important
/// class in the library, since it defines fundamental attributes of a general
/// chemical species such as its elemental formula, electrical charge and molar
/// mass. In addition, it provides the functionality to calculate its standard
/// chemical potential at given temperature *T* and pressure *P* points.
/// @see Phase, Phases
/// @ingroup Core
class Species
{
public:
    /// Construct a default Species instance
    Species();

    /// Construct a copy of a Species instance
    Species(const Species& other);

    /// Destroy this Species instance
    virtual ~Species();

    /// Assign a Species instance to this instance
    auto operator=(Species other) -> Species&;

    /// Set the name of the species
    auto setName(const std::string& name) -> Species&;

    /// Set the chemical formula of the species
    auto setFormula(const std::string& formula) -> Species&;

    /// Set the names of the elements that compose the species
    auto setElementNames(const std::vector<std::string>& element_names) -> Species&;

    /// Set the number of atoms of the elements that compose the species
    auto setElementAtoms(const std::vector<double>& element_atoms) -> Species&;

    /// Set the molar mass of the species (in units of kg/mol)
    auto setMolarMass(double val) -> Species&;

    /// Set the electrical charge of the species
    auto setCharge(double val) -> Species&;

    /// Set the thermodynamic model of the species
    auto setThermoModel(const SpeciesThermoModel& thermo_model) -> Species&;

    /// Get the name of the species
    auto name() const -> const std::string&;

    /// Get the chemical formula of the species
    auto formula() const -> const std::string&;

    /// Get the names of the elements that compose the species
    auto elementNames() const -> const std::vector<std::string>&;

    /// Get the numer of atoms of the elements that compose the species
    auto elementAtoms() const -> const std::vector<double>&;

    /// Get the molar mass of the species (in units of kg/mol)
    auto molarMass() const -> double;

    /// Get the electrical charge of the species
    auto charge() const -> double;

    /// Get the thermodynamic model of the species
    auto thermoModel() const -> const SpeciesThermoModel&;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

/// A type used to describe the thermodynamic model of a species
struct SpeciesThermoModel
{
    /// The apparent standard molar Gibbs free energy function of the species (in units of J/mol).
    ThermoPropertyFunction gibbs_energy;

    /// The apparent standard molar Helmholtz free energy function of the species (in units of J/mol).
    ThermoPropertyFunction helmholtz_energy;

    /// The apparent standard molar internal energy function of the species (in units of J/mol).
    ThermoPropertyFunction internal_energy;

    /// The apparent standard molar enthalpy function of the species (in units of J/mol).
    ThermoPropertyFunction enthalpy;

    /// The standard molar entropy function of the species (in units of J/K).
    ThermoPropertyFunction entropy;

    /// The standard molar volume function of the species (in units of m3/mol).
    ThermoPropertyFunction volume;

    /// The standard molar isobaric heat capacity function of the species (in units of J/(mol K)).
    ThermoPropertyFunction heat_capacity_cp;
};

} // namespace Reaktor
