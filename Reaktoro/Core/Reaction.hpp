// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

// C++ includes
#include <functional>
#include <memory>
#include <string>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Math/Matrix.hpp>
#include <Reaktoro/Common/ReactionEquation.hpp>
#include <Reaktoro/Common/ScalarTypes.hpp>
#include <Reaktoro/Core/Species.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalSystem;
class ChemicalProperties;

/// The function signature of the rate of a reaction (in units of mol/s).
/// @param properties The thermodynamic properties of the chemical system at (*T*, *P*, **n**)
/// @return The rate of the reaction and its partial derivatives (in units of mol/s)
/// @see Reaction
/// @ingroup Core
using ReactionRateFunction = std::function<ChemicalScalar(const ChemicalProperties&)>;

/// The function signature of the rates of a collection of reactions (in units of mol/s).
/// @param properties The thermodynamic properties of the chemical system at (*T*, *P*, **n**)
/// @see Reaction
/// @ingroup Core
using ReactionRateVectorFunction = std::function<ChemicalVector(const ChemicalProperties&)>;

/// Provide a computational representation of a chemical reaction.
/// The Reaction class provides a representation of a chemical reaction
/// and operations such as the calculation of equilibrium constants at
/// given temperature and pressure points, reaction quotients, and
/// reaction rates.
/// @see ReactionRate, EquilibriumConstant
/// @ingroup Core
class Reaction
{
public:
    /// Construct a default Reaction instance
    Reaction();

    /// Construct a Reaction instance from a ReactionEquation instance
    Reaction(const ReactionEquation& equation, const ChemicalSystem& system);

    /// Construct a copy of a Reaction instance
    Reaction(const Reaction& other);

    /// Destroy this instance
    virtual ~Reaction();

    /// Assign a Reaction instance to this instance
    auto operator=(Reaction other) -> Reaction&;

    /// Set the name of the reaction.
    auto setName(std::string name) -> void;

    /// Set the equilibrium constant function of the reaction (in natural log scale).
    auto setEquilibriumConstant(const ThermoScalarFunction& lnk) -> void;

    /// Set the rate function of the reaction (in units of mol/s).
    auto setRate(const ReactionRateFunction& function) -> void;

    /// Return the name of the reaction.
    auto name() const -> std::string;

    /// Return the equilibrium constant function of the reaction.
    auto equilibriumConstant() const -> const ThermoScalarFunction&;

    /// Return the rate function of the reaction.
    auto rate() const -> const ReactionRateFunction&;

    /// Return the equation of the reaction
    auto equation() const -> const ReactionEquation&;

    /// Return the chemical system instance of the reaction
    auto system() const -> const ChemicalSystem&;

    /// Return the reacting species of the reaction
    auto species() const -> const std::vector<Species>&;

    /// Return the indices of the reacting species of the reaction
    auto indices() const -> const Indices&;

    /// Return the stoichiometries of the reacting species of the reaction
    auto stoichiometries() const -> VectorConstRef;

    /// Return the stoichiometry of a species in the reaction equation.
    /// @param species The name of the species.
    auto stoichiometry(std::string species) const -> double;

    /// Calculate the equilibrium constant of the reaction (in natural log).
    /// @param properties The chemical properties of the system
    auto lnEquilibriumConstant(const ChemicalProperties& properties) const -> ThermoScalar;

    /// Calculate the reaction quotient of the reaction (in natural log scale).
    /// The reaction quotient of a reaction is defined as:
    /// @f[\ln Q=\sum_{i=1}^{N}\nu_{i}\ln a_{i},@f]
    /// where @f$N@f$ denotes the number of species in the multiphase system,
    /// @f$a_{i}@f$ the activity of the @f$i@f$-th species, and
    /// @f$\nu_{i}@f$ the stoichiometry of the @f$i@f$-th species in the reaction:
    /// @f[0\rightleftharpoons\sum_{i=1}^{N}\nu_{i}\alpha_{i},@f]
    /// with @f$\alpha_{i}@f$ denoting the @f$i@f$-th species. The sign
    /// convention for the stoichiometric coefficients is: *positive* for
    /// products, *negative* for reactants.
    /// @param properties The chemical properties of the system
    auto lnReactionQuotient(const ChemicalProperties& properties) const -> ChemicalScalar;

    /// Calculate the equilibrium index of the reaction as @f$\ln(Q/K)@f$.
    auto lnEquilibriumIndex(const ChemicalProperties& properties) const -> ChemicalScalar;

    /// Calculate the rate of the reaction (in units of mol/s).
    /// @param properties The thermodynamic properties of the chemical system at (*T*, *P*, **n**)
    auto rate(const ChemicalProperties& properties) const -> ChemicalScalar;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

/// Compare two Reaction instances for less than
auto operator<(const Reaction& lhs, const Reaction& rhs) -> bool;

/// Compare two Reaction instances for equality
auto operator==(const Reaction& lhs, const Reaction& rhs) -> bool;

} // namespace Reaktoro
