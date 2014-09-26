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
class  Phase;
class  Species;
struct ThermoProperties;
class  ThermoVector;
class  ThermoScalar;

/// Get the number of species in a phase
/// @param phase The phase instance
auto numSpecies(const Phase& phase) -> unsigned;

/// Get the index of a species in a phase
/// @param phase The phase instance
/// @param name The name of the species
/// @return The index of the species if found. Otherwise the number of species in the phase
auto indexSpecies(const Phase& phase, const std::string& name) -> Index;

/// Check if a phase contains a species
/// @param phase The phase instance
/// @param name The name of the species
auto containsSpecies(const Phase& phase, const std::string& name) -> bool;

/// Get the names of the phases in a container of phases
/// @param phases The container of phases
auto phaseNames(const std::vector<Phase>& phases) -> std::vector<std::string>;

/// Calculate the standard molar volumes of the species in a phase (in units of m3/mol)
/// @param phase The phase instance
/// @param T The temperature (in units of K)
/// @param P The pressure (in units of Pa)
auto volumes(const Phase& phase, double T, double P) -> ThermoProperties;

/// Calculate the standard molar entropies of the species in a phase (in units of J/K)
/// @param phase The phase instance
/// @param T The temperature (in units of K)
/// @param P The pressure (in units of Pa)
auto entropies(const Phase& phase, double T, double P) -> ThermoProperties;

/// Calculate the apparent standard molar Helmholtz free energies of the species in a phase (in units of J/mol)
/// @param phase The phase instance
/// @param T The temperature (in units of K)
/// @param P The pressure (in units of Pa)
auto helmholtzEnergies(const Phase& phase, double T, double P) -> ThermoProperties;

/// Calculate the apparent standard molar internal energies of the species in a phase (in units of J/mol)
/// @param phase The phase instance
/// @param T The temperature (in units of K)
/// @param P The pressure (in units of Pa)
auto internalEnergies(const Phase& phase, double T, double P) -> ThermoProperties;

/// Calculate the apparent standard molar enthalpies of the species in a phase (in units of J/mol)
/// @param phase The phase instance
/// @param T The temperature (in units of K)
/// @param P The pressure (in units of Pa)
auto enthalpies(const Phase& phase, double T, double P) -> ThermoProperties;

/// Calculate the apparent standard molar Gibbs free energies of the species in a phase (in units of J/mol)
/// @param phase The phase instance
/// @param T The temperature (in units of K)
/// @param P The pressure (in units of Pa)
auto gibbsEnergies(const Phase& phase, double T, double P) -> ThermoProperties;

/// Calculate the standard molar isobaric heat capacities of the species in a phase (in units of J/(mol K))
/// @param phase The phase instance
/// @param T The temperature (in units of K)
/// @param P The pressure (in units of Pa)
auto heatCapacitiesCp(const Phase& phase, double T, double P) -> ThermoProperties;

/// Calculate the molar fractions of the species
/// @param phase The phase instance
/// @param n The molar amounts of the species in the phase
auto molarFractions(const Phase& phase, const Vector& n) -> ThermoVector;

/// Calculate the concentrations of the species in a phase
/// @param phase The phase instance
/// @param n The molar amounts of the species in the phase
auto concentrations(const Phase& phase, const Vector& n) -> ThermoVector;

/// Calculate the activities of the species in a phase
/// @param phase The phase instance
/// @param T The temperature (in units of K)
/// @param P The pressure (in units of Pa)
/// @param n The molar amounts of the species in the phase
auto activities(const Phase& phase, double T, double P, const Vector& n) -> ThermoVector;

/// Calculate the density of the phase (in units of kg/m3)
/// @param phase The phase instance
/// @param T The temperature (in units of K)
/// @param P The pressure (in units of Pa)
/// @param n The molar amounts of the species in the phase
auto density(const Phase& phase, double T, double P, const Vector& n) -> ThermoScalar;

} // namespace Reaktor
