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
#include <Reaktor/Common/Matrix.hpp>
#include <Reaktor/Core/Element.hpp>

namespace Reaktor {

/// A type used to define the attributes of a Species instance
/// @see Species
/// @ingroup Core
struct SpeciesData
{
    /// The name of the chemical species
    std::string name;

    /// The chemical formula of the chemical species
    std::string formula;

    /// The elements that compose the chemical species
    std::vector<Element> elements;

    /// The number of atoms of the elements that compose the chemical species
    std::vector<double> atoms;

    /// The electrical charge of the chemical species
    double charge;

    /// The molar mass of the chemical species (in units of kg/mol)
    double molar_mass;
};

/// A type used to describe a chemical species and its attributes.
/// The Species class is used to represent a chemical species. It is an important
/// class in the library, since it defines fundamental attributes of a general
/// chemical species such as its elemental formula, electrical charge and molar mass.
/// @see Phase
/// @ingroup Core
class Species
{
public:
    /// Construct a default Species instance
    Species();

    /// Construct a custom Species instance with all its attributes
    Species(const SpeciesData& data);

    /// Construct a custom Species instance with all its attributes
    Species(std::string name, std::string formula, std::vector<Element> elements, std::vector<double> atoms, double charge, double molar_mass);

    /// Get the number of elements of the chemical species
    auto numElements() const -> unsigned;

    /// Get the name of the chemical species
    auto name() const -> const std::string&;

    /// Get the formula of the chemical species
    auto formula() const -> const std::string&;

    /// Get the names of the elements that compose the chemical species
    auto elements() const -> const std::vector<Element>&;

    /// Get the number of atoms of the elements that compose the chemical species
    auto atoms() const -> const std::vector<double>&;

    /// Get the electrical charge of the chemical species
    auto charge() const -> double;

    /// Get the molar mass of the chemical species (in units of kg/mol)
    auto molarMass() const -> double;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

/// Compare two Species instances for less than
auto operator<(const Species& lhs, const Species& rhs) -> bool;

/// Compare two Species instances for equality
auto operator==(const Species& lhs, const Species& rhs) -> bool;

/// Return the number of atoms of an element in a species
/// @param element The element instance
/// @param species The species instance
auto atoms(const Element& element, const Species& species) -> double;

/// Return the list of elements (in alphabetical order) that compose a list of species
auto collectElements(const std::vector<Species>& species) -> std::vector<Element>;

/// Return the electrical charges of all species in a list of species
auto charges(const std::vector<Species>& species) -> Vector;

/// Return the molar masses of all species in a list of species (in units of kg/mol)
auto molarMasses(const std::vector<Species>& species) -> Vector;

} // namespace Reaktor
