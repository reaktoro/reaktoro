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
#include <string>
#include <memory>

// Reaktoro includes
#include <Reaktoro/Common/ScalarTypes.hpp>

#ifdef REAKTORO_USING_THERMOFUN
// Forward declarations for ThermoFun
namespace ThermoFun { class Database; }
#endif

namespace Reaktoro {

// Forward declarations
class Database;
struct SpeciesThermoState;
struct WaterThermoState;

/// A type to calculate thermodynamic properties of chemical species
class Thermo
{
public:
    /// Construct a Thermo instance with given Database instance
    explicit Thermo(const Database& database);

#ifdef REAKTORO_USING_THERMOFUN
    /// Construct a Thermo instance with given ThermoFun::Database instance
    explicit Thermo(const ThermoFun::Database& database);
#endif

    /// Calculate the apparent standard molar Gibbs free energy of a species (in units of J/mol).
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    /// @param species The name of the species
    auto standardPartialMolarGibbsEnergy(double T, double P, std::string species) const -> ThermoScalar;

    /// Calculate the apparent standard molar Helmholtz free energy of a species (in units of J/mol).
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    /// @param species The name of the species
    auto standardPartialMolarHelmholtzEnergy(double T, double P, std::string species) const -> ThermoScalar;

    /// Calculate the apparent standard molar internal energy of a species (in units of J/mol).
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    /// @param species The name of the species
    auto standardPartialMolarInternalEnergy(double T, double P, std::string species) const -> ThermoScalar;

    /// Calculate the apparent standard molar enthalpy of a species (in units of J/mol).
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    /// @param species The name of the species
    auto standardPartialMolarEnthalpy(double T, double P, std::string species) const -> ThermoScalar;

    /// Calculate the standard molar entropies of a species (in units of J/K).
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    /// @param species The name of the species
    auto standardPartialMolarEntropy(double T, double P, std::string species) const -> ThermoScalar;

    /// Calculate the standard molar volumes of a species (in units of m3/mol).
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    /// @param species The name of the species
    auto standardPartialMolarVolume(double T, double P, std::string species) const -> ThermoScalar;

    /// Calculate the standard molar isobaric heat capacity of a species (in units of J/(mol*K)).
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    /// @param species The name of the species
    auto standardPartialMolarHeatCapacityConstP(double T, double P, std::string species) const -> ThermoScalar;

    /// Calculate the standard molar isochoric heat capacity of a species (in units of J/(mol*K)).
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    /// @param species The name of the species
    auto standardPartialMolarHeatCapacityConstV(double T, double P, std::string species) const -> ThermoScalar;

    /// Calculate the ln equilibrium constant of a reaction.
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    /// @param reaction The reaction equation
    auto lnEquilibriumConstant(double T, double P, std::string reaction) -> ThermoScalar;

    /// Calculate the log equilibrium constant of a reaction.
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    /// @param reaction The reaction equation
    auto logEquilibriumConstant(double T, double P, std::string reaction) -> ThermoScalar;

    /// Return true if there is support for the calculation of the apparent standard molar Gibbs free energy of a species.
    /// @param species The name of the species
    auto hasStandardPartialMolarGibbsEnergy(std::string species) const -> bool;

    /// Return true if there is support for the calculation of the apparent standard molar Helmholtz free energy of a species.
    /// @param species The name of the species
    auto hasStandardPartialMolarHelmholtzEnergy(std::string species) const -> bool;

    /// Return true if there is support for the calculation of the apparent standard molar internal energy of a species.
    /// @param species The name of the species
    auto hasStandardPartialMolarInternalEnergy(std::string species) const -> bool;

    /// Return true if there is support for the calculation of the apparent standard molar enthalpy of a species.
    /// @param species The name of the species
    auto hasStandardPartialMolarEnthalpy(std::string species) const -> bool;

    /// Return true if there is support for the calculation of the standard molar entropies of a species.
    /// @param species The name of the species
    auto hasStandardPartialMolarEntropy(std::string species) const -> bool;

    /// Return true if there is support for the calculation of the standard molar volumes of a species.
    /// @param species The name of the species
    auto hasStandardPartialMolarVolume(std::string species) const -> bool;

    /// Return true if there is support for the calculation of the standard molar isobaric heat capacity of a species.
    /// @param species The name of the species
    auto hasStandardPartialMolarHeatCapacityConstP(std::string species) const -> bool;

    /// Return true if there is support for the calculation of the standard molar isochoric heat capacity of a species.
    /// @param species The name of the species
    auto hasStandardPartialMolarHeatCapacityConstV(std::string species) const -> bool;

    /// Calculate the thermodynamic state of an aqueous species using the HKF model.
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    /// @param species The name of the species
    /// @see SpeciesThermoState
    auto speciesThermoStateHKF(double T, double P, std::string species) -> SpeciesThermoState;

    /// Calculate the thermodynamic state of water using the Haar--Gallagher--Kell (1984) equation of state.
    /// @param T The temperature of water (in units of K)
    /// @param P The pressure of water (in units of Pa)
    /// @see WaterThermoState
    auto waterThermoStateHGK(double T, double P) -> WaterThermoState;

    /// Calculate the thermodynamic state of water using the Wagner and Pruss (1995) equation of state.
    /// @param T The temperature of water (in units of K)
    /// @param P The pressure of water (in units of Pa)
    /// @see WaterThermoState
    auto waterThermoStateWagnerPruss(double T, double P) -> WaterThermoState;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

} // namespace Reaktoro
