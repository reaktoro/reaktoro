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
#include <vector>
#include <string>

// Reaktor includes
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Common/Matrix.hpp>
#include <Reaktor/Core/Reactions.hpp>

namespace Reaktor {

//// Forward declarations
// class Multiphase;
// class Reaction;
//struct ReactionThermoModel;
//struct ThermoScalar;
//struct ThermoVector;
// class ChemicalScalar;
// class ChemicalVector;
//
///// Return the number of species in a reaction
///// @param reaction The reaction instance
//auto numSpecies(const Reaction& reaction) -> unsigned;
//
///// Check if a reaction contains a species
///// @param reaction The reaction instance
///// @param species The name of the species
//auto containsSpecies(const Reaction& reaction, const std::string& species) -> bool;
//
///// Return the index of a species in a reaction
///// @param reaction The reaction instance
///// @param species The name of the species
//auto indexSpecies(const Reaction& reaction, const std::string& species) -> Index;
//
///// Return the indices of the phases that participates in a reaction
///// @param multiphase The multiphase system
///// @param reaction The reaction
//auto phaseIndicesInReaction(const Multiphase& multiphase, const Reaction& reaction) -> Indices;
//
///// Return the indices of the reactions that contains a species
///// @param reactions The set of reactions
///// @param ispecies The index of the species
//auto indicesReactionsWithSpecies(const Reactions& reactions, const Index& ispecies) -> Indices;
//
///// Return the stoichiometry of a species in a reaction
///// @param reaction The reaction instance
///// @param species The name of the species
///// @return The stoichiometry of the species if it participates
///// in the reaction, or zero otherwise.
//auto stoichiometry(const Reaction& reaction, const std::string& species) -> double;
//
///// Create a ReactionThermoModel instance for a reaction.
///// The created thermodynamic model for the reaction uses the thermodynamic models of
///// its reacting species.
///// @param multiphase The multiphase instance
///// @param reaction The reaction instance
//auto thermoModel(const Multiphase& multiphase, const Reaction& reaction) -> ReactionThermoModel;
//
///// Calculate the natural logarithm of the equilibrium constant of a reaction
///// @param reaction The reaction instance
///// @param T The temperature of the chemical system (in units of K)
///// @param P The pressure of the chemical system (in units of Pa)
//auto equilibriumConstant(const Reaction& reaction, double T, double P) -> ThermoScalar;
//
///// Calculate the natural logarithm of the equilibrium constants of a set of reactions
///// @param reactions The set of reactions
///// @param T The temperature of the chemical system (in units of K)
///// @param P The pressure of the chemical system (in units of Pa)
//auto equilibriumConstants(const Reactions& reactions, double T, double P) -> ThermoVector;
//
///// Calculate the standard molar entropy of a reaction (in units of J/K)
///// @param species The species instance
///// @param T The temperature for the calculation (in units of K)
///// @param P The pressure for the calculation (in units of Pa)
//auto entropy(const Reaction& reaction, double T, double P) -> ThermoScalar;
//
///// Calculate the standard molar entropy of a set of reactions (in units of J/K)
///// @param reactions The set of reactions
///// @param T The temperature for the calculation (in units of K)
///// @param P The pressure for the calculation (in units of Pa)
//auto entropies(const Reactions& reactions, double T, double P) -> ThermoVector;
//
///// Calculate the apparent standard molar Helmholtz free energy of a reaction (in units of J/mol)
///// @param species The species instance
///// @param T The temperature for the calculation (in units of K)
///// @param P The pressure for the calculation (in units of Pa)
//auto helmholtzEnergy(const Reaction& reaction, double T, double P) -> ThermoScalar;
//
///// Calculate the apparent standard molar Helmholtz free energy of a set of reactions (in units of J/mol)
///// @param reactions The set of reactions
///// @param T The temperature for the calculation (in units of K)
///// @param P The pressure for the calculation (in units of Pa)
//auto helmholtzEnergies(const Reactions& reactions, double T, double P) -> ThermoVector;
//
///// Calculate the apparent standard molar internal energy of a reaction (in units of J/mol)
///// @param species The species instance
///// @param T The temperature for the calculation (in units of K)
///// @param P The pressure for the calculation (in units of Pa)
//auto internalEnergy(const Reaction& reaction, double T, double P) -> ThermoScalar;
//
///// Calculate the apparent standard molar internal energy of a set of reactions (in units of J/mol)
///// @param reactions The set of reactions
///// @param T The temperature for the calculation (in units of K)
///// @param P The pressure for the calculation (in units of Pa)
//auto internalEnergies(const Reactions& reactions, double T, double P) -> ThermoVector;
//
///// Calculate the apparent standard molar enthalpy of a reaction (in units of J/mol)
///// @param species The species instance
///// @param T The temperature for the calculation (in units of K)
///// @param P The pressure for the calculation (in units of Pa)
//auto enthalpy(const Reaction& reaction, double T, double P) -> ThermoScalar;
//
///// Calculate the apparent standard molar enthalpy of a set of reactions (in units of J/mol)
///// @param reactions The set of reactions
///// @param T The temperature for the calculation (in units of K)
///// @param P The pressure for the calculation (in units of Pa)
//auto enthalpies(const Reactions& reactions, double T, double P) -> ThermoVector;
//
///// Calculate the apparent standard molar Gibbs free energy of a reaction (in units of J/mol)
///// @param species The species instance
///// @param T The temperature for the calculation (in units of K)
///// @param P The pressure for the calculation (in units of Pa)
//auto gibbsEnergy(const Reaction& reaction, double T, double P) -> ThermoScalar;
//
///// Calculate the apparent standard molar Gibbs free energy of a set of reactions (in units of J/mol)
///// @param reactions The set of reactions
///// @param T The temperature for the calculation (in units of K)
///// @param P The pressure for the calculation (in units of Pa)
//auto gibbsEnergies(const Reactions& reactions, double T, double P) -> ThermoVector;
//
///// Calculate the kinetic rate of a reaction
///// @param reaction The reaction instance
///// @param T The temperature of the chemical system (in units of K)
///// @param P The pressure of the chemical system (in units of Pa)
///// @param n The molar abundance of the species in the chemical system (in units of mol)
///// @param a The activities of every species in the chemical system and their molar derivatives
///// @return The rate of the reaction and its molar derivatives
//auto rate(const Reaction& reaction, double T, double P, const Vector& n, const ChemicalVector& a) -> ChemicalScalar;
//
///// Calculate the kinetic rates of a set of reactions
///// @param reactions The set of reactions
///// @param T The temperature of the chemical system (in units of K)
///// @param P The pressure of the chemical system (in units of Pa)
///// @param n The molar abundance of the species in the chemical system (in units of mol)
///// @param a The activities of every species in the chemical system and their molar derivatives
///// @return The rate of the reaction and its molar derivatives
//auto rates(const Reactions& reactions, double T, double P, const Vector& n, const ChemicalVector& a) -> ChemicalVector;
//
///// Calculate the reaction quotient of the reaction.
///// The reaction quotient @f$ Q @f$ of a reaction is defined as:
///// @f[
/////     Q=\prod_{i=1}^{N}a_{i}^{\nu_{i}},
///// @f]
///// where @f$ N @f$ denotes the number of species in the chemical system,
///// @f$ a_{i} @f$ the activity of the @f$ i @f$-th species, and
///// @f$ \nu_{i} @f$ the stoichiometry of the @f$ i @f$-th species in the
///// reaction:
///// @f[
/////     0\rightleftharpoons\sum_{i=1}^{N}\nu_{i}\alpha_{i},
///// @f]
///// with @f$ \alpha_{i} @f$ denoting the @f$ i @f$-th species. The sign
///// convention for the stoichiometric coefficients is: *positive* for
///// products, *negative* for reactants.
///// @param reaction The reaction instance
///// @param a The activities of every species in the chemical system and their molar derivatives
//auto reactionQuotient(const Reaction& reaction, const ChemicalVector& a) -> ChemicalScalar;
//
///// Calculate the reaction quotients of a set of reactions
//auto reactionQuotients(const Reactions& reactions, const ChemicalVector& a) -> ChemicalVector;
//
///// Assemble the stoichiometric matrix of a set of reactions in a multiphase system
///// @param multiphase The multiphase system
///// @param reaction The set of reactions
//auto stoichiometricMatrix(const Multiphase& multiphase, const Reactions& reactions) -> Matrix;

} // namespace Reaktor
