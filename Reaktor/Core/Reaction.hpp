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
#include <functional>
#include <memory>
#include <string>
#include <vector>

// Reaktor includes
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Common/Matrix.hpp>
#include <Reaktor/Common/ReactionEquation.hpp>
#include <Reaktor/Common/ThermoScalar.hpp>
#include <Reaktor/Common/ChemicalScalar.hpp>
#include <Reaktor/Core/Species.hpp>

namespace Reaktor {

// Forward declarations
class ChemicalScalar;
class ChemicalSystem;
class ChemicalVector;

/// The function signature of the rate of a reaction (in units of mol/s).
/// @param T The temperature value (in units of K)
/// @param P The pressure value (in units of Pa)
/// @param n The molar amounts of all species in the system (in units of mol)
/// @param ln_a The activities of all species in the system and their molar derivatives (in natural log scale)
/// @return The rate of the reaction and its molar derivatives (in units of mol/s)
/// @see Reaction
/// @ingroup Core
typedef std::function<
    ChemicalScalar(double, double, const Vector&, const ChemicalVector&)>
        ReactionRateFunction;

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

    /// Return a copy of this Reaction instance with a new equilibrium constant function
    auto withEquilibriumConstant(const ThermoScalarFunction& lnk) const -> Reaction;

    /// Return a copy of this Reaction instance with a new reaction rate function
    auto withRate(const ChemicalScalarFunction& rate) const -> Reaction;

    /// Return the equation of the reaction
    auto equation() const -> const ReactionEquation&;

    /// Return the chemical system instance of the reaction
    auto system() const -> const ChemicalSystem&;

    /// Return the reacting species of the reaction
    auto species() const -> const std::vector<Species>&;

    /// Return the indices of the reacting species of the reaction
    auto indices() const -> const Indices&;

    /// Return the stoichiometries of the reacting species of the reaction
    auto stoichiometries() const -> const Vector&;

    /// Calculate the equilibrium constant of the reaction (in natural log).
    auto lnEquilibriumConstant(double T, double P) const -> ThermoScalar;

    /// Calculate the apparent standard molar Gibbs free energy of the reaction (in units of J/mol).
    auto standardGibbsEnergy(double T, double P) const -> ThermoScalar;

    /// Calculate the apparent standard molar Helmholtz free energy of the reaction (in units of J/mol).
    auto standardHelmholtzEnergy(double T, double P) const -> ThermoScalar;

    /// Calculate the apparent standard molar internal energy of the reaction (in units of J/mol).
    auto standardInternalEnergy(double T, double P) const -> ThermoScalar;

    /// Calculate the apparent standard molar enthalpy of the reaction (in units of J/mol).
    auto standardEnthalpy(double T, double P) const -> ThermoScalar;

    /// Calculate the standard molar entropies of the reaction (in units of J/K).
    auto standardEntropy(double T, double P) const -> ThermoScalar;

    /// Calculate the standard molar volumes of the reaction (in units of m3/mol).
    auto standardVolume(double T, double P) const -> ThermoScalar;

    /// Calculate the standard molar isobaric heat capacity of the reaction (in units of J/(mol*K)).
    auto standardHeatCapacity(double T, double P) const -> ThermoScalar;

    /// Calculate the rate of the reaction (in units of mol/s).
    auto rate(double T, double P, const Vector& n) const -> ChemicalScalar;

    /// Calculate the reaction quotient of the reaction (in natural log scale).
    /// The reaction quotient of a reaction is defined as:
    /// @f[
    ///     \ln Q_r=\sum_{i=1}^{N}\nu_{i}\ln a_{i},
    /// @f]
    /// where @f$N@f$ denotes the number of species in the chemical system,
    /// @f$a_{i}@f$ the activity of the @f$i@f$-th species, and
    /// @f$\nu_{i}@f$ the stoichiometry of the @f$i@f$-th species in the
    /// reaction:
    /// @f[
    ///     0\rightleftharpoons\sum_{i=1}^{N}\nu_{i}\alpha_{i},
    /// @f]
    /// with @f$\alpha_{i}@f$ denoting the @f$i@f$-th species. The sign
    /// convention for the stoichiometric coefficients is: *positive* for
    /// products, *negative* for reactants.
    auto lnReactionQuotient(double T, double P, const Vector& n) const -> ChemicalScalar;

    /// Calculate the equilibrium index of the reaction (in natural log scale).
    /// The equilibrium index of a reaction is defined as @f$\ln\Omega_r=\ln Q_r-\ln K_r@f$,
    /// where @f$\ln Q_r@f$ is the reaction quotient and @f$\ln K_r@f$ is the equilibrium
    /// constant of the reaction. For mineral reactions, the equilibrium index is exactly its
    /// saturation index.
    auto lnEquilibriumIndex(double T, double P, const Vector& n) const -> ChemicalScalar;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

} // namespace Reaktor
