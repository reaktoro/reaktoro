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
#include <string>
#include <vector>

// Reaktor includes
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Common/Vector.hpp>

namespace Reaktor {

// Forward declarations
class  Species;
struct ThermoScalar;
struct ThermoVector;

///// Get the number of elements in a species
///// @param species The species instance
//auto numElements(const Species& species) -> unsigned;
//
///// Check if a species contains a chemical element
///// @param species The species instance
///// @param element The name of the chemical element
//auto containsElement(const Species& species, const std::string& element) -> bool;
//
///// Get the index of a chemical element in a species
///// @param species The species instance
///// @param element The name of the chemical element
///// @param The index of the chemical element if found, or the number of elements otherwise
//auto elementIndex(const Species& species, const std::string& element) -> Index;
//
///// Get the number of atoms of an element in the species
///// @param species The species instance
///// @param element The name of the chemical element
//auto elementAtoms(const Species& species, const std::string& element) -> double;
//
///// Calculate the standard molar volume of a species (in units of m3/mol)
///// @param species The species instance
///// @param T The temperature (in units of K)
///// @param P The pressure (in units of Pa)
//auto volume(const Species& species, double T, double P) -> ThermoScalar;
//
///// Calculate the standard molar volumes of a collection of species (in units of m3/mol)
///// @param species The collection of species
///// @param T The temperature (in units of K)
///// @param P The pressure (in units of Pa)
//auto volumes(const std::vector<Species>& species, double T, double P) -> ThermoVector;
//
///// Calculate the standard molar entropy of the species (in units of J/K)
///// @param species The species instance
///// @param T The temperature (in units of K)
///// @param P The pressure (in units of Pa)
//auto entropy(const Species& species, double T, double P) -> ThermoScalar;
//
///// Calculate the standard molar entropy of the species (in units of J/K)
///// @param species The collection of species
///// @param T The temperature (in units of K)
///// @param P The pressure (in units of Pa)
//auto entropies(const std::vector<Species>& species, double T, double P) -> ThermoVector;
//
///// Calculate the apparent standard molar Helmholtz free energy of the species (in units of J/mol)
///// @param species The species instance
///// @param T The temperature (in units of K)
///// @param P The pressure (in units of Pa)
//auto helmholtzEnergy(const Species& species, double T, double P) -> ThermoScalar;
//
///// Calculate the apparent standard molar Helmholtz free energy of the species (in units of J/mol)
///// @param species The collection of species
///// @param T The temperature (in units of K)
///// @param P The pressure (in units of Pa)
//auto helmholtzEnergies(const std::vector<Species>& species, double T, double P) -> ThermoVector;
//
///// Calculate the apparent standard molar internal energy of the species (in units of J/mol)
///// @param species The species instance
///// @param T The temperature (in units of K)
///// @param P The pressure (in units of Pa)
//auto internalEnergy(const Species& species, double T, double P) -> ThermoScalar;
//
///// Calculate the apparent standard molar internal energy of the species (in units of J/mol)
///// @param species The collection of species
///// @param T The temperature (in units of K)
///// @param P The pressure (in units of Pa)
//auto internalEnergies(const std::vector<Species>& species, double T, double P) -> ThermoVector;
//
///// Calculate the apparent standard molar enthalpy of the species (in units of J/mol)
///// @param species The species instance
///// @param T The temperature (in units of K)
///// @param P The pressure (in units of Pa)
//auto enthalpy(const Species& species, double T, double P) -> ThermoScalar;
//
///// Calculate the apparent standard molar enthalpy of the species (in units of J/mol)
///// @param species The collection of species
///// @param T The temperature (in units of K)
///// @param P The pressure (in units of Pa)
//auto enthalpies(const std::vector<Species>& species, double T, double P) -> ThermoVector;
//
///// Calculate the apparent standard molar Gibbs free energy of the species (in units of J/mol)
///// @param species The species instance
///// @param T The temperature (in units of K)
///// @param P The pressure (in units of Pa)
//auto gibbsEnergy(const Species& species, double T, double P) -> ThermoScalar;
//
///// Calculate the apparent standard molar Gibbs free energy of the species (in units of J/mol)
///// @param species The collection of species
///// @param T The temperature (in units of K)
///// @param P The pressure (in units of Pa)
//auto gibbsEnergies(const std::vector<Species>& species, double T, double P) -> ThermoVector;
//
///// Calculate the standard molar isobaric heat capacity of the species (in units of J/(mol K))
///// @param species The species instance
///// @param T The temperature (in units of K)
///// @param P The pressure (in units of Pa)
//auto heatCapacityCp(const Species& species, double T, double P) -> ThermoScalar;
//
///// Calculate the standard molar isobaric heat capacity of the species (in units of J/(mol K))
///// @param species The collection of species
///// @param T The temperature (in units of K)
///// @param P The pressure (in units of Pa)
//auto heatCapacitiesCp(const std::vector<Species>& species, double T, double P) -> ThermoVector;
//
///// Get the names of the species in a container of species
///// @param The container of species
//auto speciesNames(const std::vector<Species>& species) -> std::vector<std::string>;
//
///// Get the electrical charges of the species in a container of species
///// @param The container of species
//auto speciesCharges(const std::vector<Species>& species) -> Vector;
//
///// Get the molar masses of the species in a container of species
///// @param The container of species
//auto speciesMolarMasses(const std::vector<Species>& species) -> Vector;

} // namespace Reaktor
