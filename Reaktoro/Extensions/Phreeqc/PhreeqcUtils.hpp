// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/AggregateState.hpp>

// Phreeqc includes
#include <Reaktoro/Extensions/Phreeqc/PhreeqcLegacy.hpp>

//==================================================
// WARNING WARNING WARNING WARNING WARNING WARNING
//==================================================
// This header file must not be included by another
// header file that will be exposed to users.
// If so, this propagates the need for phreeqc
// header files to available in the user system.
//==================================================

namespace Reaktoro {
namespace PhreeqcUtils {

/// Load the PHREEQC instance with a database file.
/// @param phreeqc The PHREEQC instance
/// @param database The path to the database file, including its file name, or a
/// multi-line string containing the database contents itself
auto load(PHREEQC& phreeqc, String database) -> void;

/// Execute a PHREEQC input script.
/// @param phreeqc The PHREEQC instance
/// @param input The path to the input script, including its file name, or a
/// multi-line string containing the input script contents itself
/// @param output The file name where the result is output
auto execute(PHREEQC& phreeqc, String input, String output) -> void;

/// Find an element in a Phreeqc instance.
/// @param phreeqc The Phreeqc instance
/// @param name The name of the element
/// @return A pointer to the element if found, nullptr otherwise.
auto findElement(const PHREEQC& phreeqc, String name) -> PhreeqcElement*;

/// Find a species in a Phreeqc instance.
/// @param phreeqc The Phreeqc instance
/// @param name The name of the species
/// @return A pointer to the species if found, nullptr otherwise.
auto findSpecies(const PHREEQC& phreeqc, String name) -> PhreeqcSpecies*;

/// Find a phase in a Phreeqc instance.
/// @param phreeqc The Phreeqc instance
/// @param name The name of the phase
/// @return A pointer to the phase if found, nullptr otherwise.
auto findPhase(const PHREEQC& phreeqc, String name) -> PhreeqcPhase*;

/// Return true if the Phreeqc species instance is an aqueous species.
auto isAqueousSpecies(const PhreeqcSpecies* species) -> bool;

/// Return true if the Phreeqc phase instance is a gaseous species.
auto isGaseousSpecies(const PhreeqcPhase* phase) -> bool;

/// Return true if the Phreeqc phase instance is a mineral species.
auto isMineralSpecies(const PhreeqcPhase* phase) -> bool;

/// Return true if the Phreeqc species instance is an exchange species.
auto isExchangeSpecies(const PhreeqcSpecies* species) -> bool;

/// Return true if the Phreeqc species instance is a surface species.
auto isSurfaceSpecies(const PhreeqcSpecies* species) -> bool;

/// Return the symbol of a Phreeqc element.
/// @param element The pointer to the Phreeqc element
auto symbol(const PhreeqcElement* species) -> String;

/// Return the name of a Phreeqc element.
/// @param element The pointer to the Phreeqc element
auto name(const PhreeqcElement* species) -> String;

/// Return the molar mass of a Phreeqc element (in kg/mol).
/// @param element The pointer to the Phreeqc element
auto molarMass(const PhreeqcElement* species) -> double;

/// Return the name of a Phreeqc phase.
/// @param phase The pointer to the Phreeqc phase
auto name(const PhreeqcPhase* phase) -> String;

/// Return the name of a Phreeqc species.
/// @param species The pointer to the Phreeqc species
auto name(const PhreeqcSpecies* species) -> String;

/// Return the name of a Phreeqc phase.
/// @param phase The pointer to the Phreeqc phase
auto name(const PhreeqcPhase* phase) -> String;

/// Return the formula of a Phreeqc species.
/// @param species The pointer to the Phreeqc species
auto formula(const PhreeqcSpecies* species) -> String;

/// Return the formula of a Phreeqc phase.
/// @param phase The pointer to the Phreeqc phase
auto formula(const PhreeqcPhase* phase) -> String;

/// Return the element symbols and their coefficients in a Phreeqc species.
/// @param species The pointer to the Phreeqc species
auto elements(const PhreeqcSpecies* species) -> Map<PhreeqcElement*, double>;

/// Return the element symbols and their coefficients in a Phreeqc phase.
/// @param phase The pointer to the Phreeqc phase
auto elements(const PhreeqcPhase* phase) -> Map<PhreeqcElement*, double>;

/// Return the charge of a Phreeqc species.
/// @param species The pointer to the Phreeqc species
auto charge(const PhreeqcSpecies* species) -> double;

/// Return the charge of a Phreeqc phase.
/// @param phase The pointer to the Phreeqc phase
auto charge(const PhreeqcPhase* phase) -> double;

/// Return the aggregate state of a Phreeqc species.
/// @param species The pointer to the Phreeqc species
auto aggregateState(const PhreeqcSpecies* species) -> AggregateState;

/// Return the aggregate state of a Phreeqc phase.
/// @param phase The pointer to the Phreeqc phase
auto aggregateState(const PhreeqcPhase* phase) -> AggregateState;

/// Return the stoichiometric coefficient of a species with given name in a PHREEQC species.
/// @param name The name of the species.
/// @param species The pointer to the Phreeqc species
auto stoichiometry(const PhreeqcSpecies* species, String name) -> double;

/// Return the stoichiometric coefficient of a species with given name in a PHREEQC phase.
/// @param name The name of the species.
/// @param species The pointer to the Phreeqc phase
auto stoichiometry(const PhreeqcPhase* phase, String name) -> double;

/// Return the reaction equation of a Phreeqc species.
/// The equation is defined by pairs of species names and their stoichiometric coefficients.
/// @note If a primary species, an empty equation is returned.
/// @param species The pointer to the Phreeqc species
auto reactionEquation(const PhreeqcSpecies* species) -> Pairs<String, double>;

/// Return the reaction equation of a Phreeqc phase.
/// The equation is defined by pairs of species names and their stoichiometric coefficients.
/// @param phase The pointer to the Phreeqc phase
auto reactionEquation(const PhreeqcPhase* phase) -> Pairs<String, double>;

/// Return the reactants in a PHREEQC species.
/// For example, for the aqueous species `CH4`, with
/// reaction `CO3-2 + 10 H+ + 8 e- = CH4 + 3 H2O`, the returned
/// pairs are: `{ {"CO3-2", 1}, {"H+", 10}, {"e-", 8}, {"H2O", -3} }`.
/// @param species The pointer to the Phreeqc species
auto reactants(const PhreeqcSpecies* species) -> Pairs<PhreeqcSpecies*, double>;

/// Return the reactants and their stoichiometric coefficients in a PHREEQC phase.
/// For example, for the mineral species `K-feldspar`, with
/// reaction `KAlSi3O8 + 8 H2O = K+ + Al(OH)4- + 3 H4SiO4`, the returned
/// pairs are: `{ {"K+", 1}, {"Al(OH)4-", 1}, {"H4SiO4", 3}, {"H2O", -8} }`.
/// @param phase The pointer to the Phreeqc phase
auto reactants(const PhreeqcPhase* phase) -> Pairs<PhreeqcSpecies*, double>;

/// Return the index of a Phreeqc species in a set of species.
/// @param name The name of the Phreeqc species
/// @param species The pointers to Phreeqc species instances
auto index(String name, const Vec<PhreeqcSpecies*>& species) -> std::size_t;

/// Return the index of a Phreeqc phase in a set of phases.
/// @param name The name of the Phreeqc phase
/// @param phases The pointers to Phreeqc phase instances
auto index(String name, const Vec<PhreeqcPhase*>& phases) -> std::size_t;

/// Return the active aqueous species in a Phreeqc instance.
/// @param phreeqc The Phreeqc instance
auto activeAqueousSpecies(const PHREEQC& phreeqc) -> Vec<PhreeqcSpecies*>;

/// Return the active exchange species in a PHREEQC instance.
auto activeExchangeSpecies(const PHREEQC& phreeqc) -> Vec<PhreeqcSpecies*>;

/// Return the active aqueous product species in a Phreeqc instance.
/// A product aqueous species is an aqueous species defined in terms of master species.
/// @param phreeqc The Phreeqc instance
auto activeProductSpecies(const PHREEQC& phreeqc) -> Vec<PhreeqcSpecies*>;

/// Return the active gaseous species in a Phreeqc instance defined in a `GAS_PHASE` block.
/// @param phreeqc The Phreeqc instance
auto activeGaseousSpecies(const PHREEQC& phreeqc) -> Vec<PhreeqcPhase*>;

/// Return the active phases in a Phreeqc instance defined in a `EQUILIBRIUM_PHASES` block.
/// @param phreeqc The Phreeqc instance
auto activePhasesInEquilibriumPhases(const PHREEQC& phreeqc) -> Vec<PhreeqcPhase*>;

/// Return the phases in a Phreeqc instance that have been marked for saturation index calculation.
/// @param phreeqc The Phreeqc instance
auto activePhasesInSaturationList(const PHREEQC& phreeqc) -> Vec<PhreeqcPhase*>;

/// Return the amounts of Phreeqc species.
/// @param phreeqc The Phreeqc instance
/// @param species The pointers to Phreeqc species instances
auto speciesAmounts(const PHREEQC& phreeqc, const Vec<PhreeqcSpecies*>& species) -> ArrayXr;

/// Return the amounts of Phreeqc phases.
/// @param phreeqc The Phreeqc instance
/// @param phases The pointers to Phreeqc phase instances
auto speciesAmounts(const PHREEQC& phreeqc, const Vec<PhreeqcPhase*>& phases) -> ArrayXr;

/// Return the equilibrium constant function (log base 10) of a Phreeqc species.
/// @param species The pointer to the Phreeqc species
auto lgEquilibriumConstantFn(const PhreeqcSpecies* species) -> Fn<real(real,real)>;

/// Return the equilibrium constant function (log base 10) of a Phreeqc phase.
/// @param phase The pointer to the Phreeqc phase
auto lgEquilibriumConstantFn(const PhreeqcPhase* phase) -> Fn<real(real,real)>;

/// Return the enthalpy of reaction function (in J/mol) of a Phreeqc species.
/// @param species The pointer to the Phreeqc species
auto enthalpyChangeFn(const PhreeqcSpecies* species) -> Fn<real(real,real)>;

/// Return the enthalpy of reaction function (in J/mol) of a Phreeqc phase.
/// @param phase The pointer to the Phreeqc phase
auto enthalpyChangeFn(const PhreeqcPhase* phase) -> Fn<real(real,real)>;

} // namespace PhreeqcUtils
} // namespace Reaktoro
