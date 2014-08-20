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
#include <iostream>
#include <string>
#include <vector>
#include <memory>

// Reaktor includes
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Common/Vector.hpp>
#include <Reaktor/Core/Functions.hpp>

namespace Reaktor {

/// Provide a computational representation of a chemical species
///
/// The Species class is used to represent a chemical species. It is an important
/// class in the library, since it defines fundamental attributes of a general
/// chemical species such as its elemental formula, electrical charge and molar
/// mass. In addition, it provides the functionality to calculate its standard
/// chemical potential at given temperature *T* and pressure *P* points.
///
/// @see Phase, ChemicalSystem
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

    /// Assigns a Species instance to this instance
    auto operator=(Species other) -> Species&;

    /// Set the name of the chemical species
    /// @return A reference to this Species instance
    auto setName(const std::string& name) -> Species&;

    /// Set the formula of the chemical species
    /// @return A reference to this Species instance
    auto setFormula(const std::string& formula) -> Species&;

    /// Set the elements that compose the chemical species
    /// @return A reference to this Species instance
    auto setElements(const std::vector<std::string>& elements) -> Species&;

    /// Set the coefficients of the elements that compose the chemical species
    /// @return A reference to this Species instance
    auto setCoefficients(const std::vector<double>& coefficients) -> Species&;

    /// Set the molar mass of the chemical species (in units of kg/mol)
    /// @return A reference to this Species instance
    auto setMolarMass(double value) -> Species&;

    /// Set the electrical charge of the chemical species
    /// @return A reference to this Species instance
    auto setCharge(double value) -> Species&;

    /// Set the standard chemical potential function of the chemical species (in units of J/mol)
    /// @return A reference to this Species instance
    auto setChemicalPotential(const ChemicalPotential& function) -> Species&;

    /// Get the name of the chemical species
    auto name() const -> const std::string&;

    /// Get the formula of the chemical species
    auto formula() const -> const std::string&;

    /// Get the elements that compose the chemical species
    auto elements() const -> const std::vector<std::string>&;

    /// Get the coefficients of the elements that compose the chemical species
    auto coefficients() const -> const std::vector<double>&;

    /// Get the molar mass of the chemical species
    auto molarMass() const -> double;

    /// Get the electrical charge of the chemical species
    auto charge() const -> double;

    /// Get the standard chemical potential function of the chemical species
    auto chemicalPotential() const -> const ChemicalPotential&;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

/// Output a Species instance
auto operator<<(std::ostream& out, const Species& species) -> std::ostream&;

} /* namespace Reaktor */
