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
#include <Reaktor/Core/Component.hpp>

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

    /// The components that compose the chemical species
    ComponentList components;

    /// The stoichiometries of the components that compose the chemical species
    std::vector<double> stoichiometries;

    /// The electrical charge of the chemical species
    double charge = 0.0;

    /// The molar mass of the chemical species (in units of kg/mol)
    double molar_mass = 0.0;
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

    /// Construct a Species instance with all its attributes
    Species(const SpeciesData& data);

    /// Get the name of the chemical species
    auto name() const -> const std::string&;

    /// Get the formula of the chemical species
    auto formula() const -> const std::string&;

    /// Get the names of the components that compose the chemical species
    auto components() const -> const ComponentList&;

    /// Get the stoichiometries of the components that compose the chemical species
    auto stoichiometries() const -> const std::vector<double>&;

    /// Get the electrical charge of the chemical species
    auto charge() const -> double;

    /// Get the molar mass of the chemical species (in units of kg/mol)
    auto molarMass() const -> double;

private:
    /// The immutable shared data of the Species class
    std::shared_ptr<SpeciesData> data;
};

/// A type used to define a list of Species instances
typedef std::vector<Species> SpeciesList;

/// Return a list of unique components that compose a collection of species
auto components(const SpeciesList& species) -> ComponentList;

/// Return the electrical charges of all species in a list of species
auto charges(const SpeciesList& species) -> std::vector<double>;

/// Return the molar masses of all species in a list of species (in units of kg/mol)
auto molarMasses(const SpeciesList& species) -> std::vector<double>;

/// Return the stoichiometry of a components in the species
/// @param component The component instance
/// @param species The species instance
auto stoichiometry(const Component& component, const Species& species) -> double;

} // namespace Reaktor
