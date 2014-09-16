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
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Common/ThermoProperty.hpp>
#include <Reaktor/Core/Functions.hpp>

namespace Reaktor {

// Forward declarations
class ThermoVector;

/// Provide a computational representation of a chemical reaction
///
/// The Reaction class provides a representation of a chemical reaction
/// and operations such as the calculation of equilibrium constants at
/// given temperature and pressure points, reaction quotients, and
/// reaction rates.
///
/// @see ReactionRate, EquilibriumConstant
/// @ingroup Core
class Reaction
{
public:
    /// Construct a default Reaction instance
    Reaction();

    /// Construct a copy of a Reaction instance
    Reaction(const Reaction& other);

    /// Destroy this instance
    virtual ~Reaction();

    /// Assign a Reaction instance to this instance
    auto operator=(Reaction other) -> Reaction&;

    /// Set the names of the reacting species of the reation
    auto setSpecies(const std::vector<std::string>& species) -> Reaction&;

    /// Set the indices of the reacting species of the reation
    auto setIndices(const Indices& indices) -> Reaction&;

    /// Set the stoichiometries of the reacting species of the reation
    auto setStoichiometries(const std::vector<double>& stoichiometries) -> Reaction&;

    /// Set the equilibrium constant function of the reaction
    auto setEquilibriumConstant(const EquilibriumConstant& eqconstant) -> Reaction&;

    /// Set the reaction rate function of the reaction
    auto setRate(const Rate& rate) -> Reaction&;

    /// Get the names of the reacting species of the reaction
    auto species() const -> const std::vector<std::string>&;

    /// Get the indices of the reacting species of the reaction
    auto indices() const -> const Indices&;

    /// Get the stoichiometries of the reacting species of the reaction
    auto stoichiometries() const -> const std::vector<double>&;

    /// Get the equilibrium constant function of the reaction
    auto equilibriumConstant() const -> const EquilibriumConstant&;

    /// Get the reaction rate function of the reaction
    auto rate() const -> const Rate&;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

struct ReactionThermoModel
{
    /// The function for the equilibrium constant of the reaction (in terms of its natural logarithm)
    ThermoPropertyFunction lnk;

    /// The function for the standard molar Gibbs free energy of the reaction  (in units of J/mol).
    ThermoPropertyFunction gibbs_energy;

    /// The function for the standard molar Helmholtz free energy of the reaction  (in units of J/mol).
    ThermoPropertyFunction helmholtz_energy;

    /// The function for the standard molar internal energy of the reaction  (in units of J/mol).
    ThermoPropertyFunction internal_energy;

    /// The function for the standard molar enthalpy of the reaction  (in units of J/mol).
    ThermoPropertyFunction enthalpy;

    /// The function for the standard molar entropy of the reaction (in units of J/K).
    ThermoPropertyFunction entropy;
};

/// Outputs the Reaction instance
auto operator<<(std::ostream& out, const Reaction& reaction) -> std::ostream&;

} // namespace Reaktor
