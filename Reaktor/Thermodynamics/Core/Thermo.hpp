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
#include <memory>

namespace Reaktor {

// Forward declarations
class Database;
class ThermoScalar;
struct SpeciesThermoState;
struct WaterThermoState;

/// A type to calculate thermodynamic properties of chemical species
class Thermo
{
public:
    /// Construct a Thermo instance with given Database instance
    explicit Thermo(const Database& database);

    /// Calculate the apparent standard molar Gibbs free energy of a species (in units of J/mol).
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    /// @param species The name of the species
    auto standardGibbsEnergy(double T, double P, std::string species) const -> ThermoScalar;

    /// Calculate the apparent standard molar Helmholtz free energy of a species (in units of J/mol).
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    /// @param species The name of the species
    auto standardHelmholtzEnergy(double T, double P, std::string species) const -> ThermoScalar;

    /// Calculate the apparent standard molar internal energy of a species (in units of J/mol).
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    /// @param species The name of the species
    auto standardInternalEnergy(double T, double P, std::string species) const -> ThermoScalar;

    /// Calculate the apparent standard molar enthalpy of a species (in units of J/mol).
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    /// @param species The name of the species
    auto standardEnthalpy(double T, double P, std::string species) const -> ThermoScalar;

    /// Calculate the standard molar entropies of a species (in units of J/K).
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    /// @param species The name of the species
    auto standardEntropy(double T, double P, std::string species) const -> ThermoScalar;

    /// Calculate the standard molar volumes of a species (in units of m3/mol).
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    /// @param species The name of the species
    auto standardVolume(double T, double P, std::string species) const -> ThermoScalar;

    /// Calculate the standard molar isobaric heat capacity of a species (in units of J/(mol*K)).
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    /// @param species The name of the species
    auto standardHeatCapacity(double T, double P, std::string species) const -> ThermoScalar;

    /// Return true if there is support for the calculation of the apparent standard molar Gibbs free energy of a species.
    /// @param species The name of the species
    auto checkStandardGibbsEnergy(std::string species) const -> bool;

    /// Return true if there is support for the calculation of the apparent standard molar Helmholtz free energy of a species.
    /// @param species The name of the species
    auto checkStandardHelmholtzEnergy(std::string species) const -> bool;

    /// Return true if there is support for the calculation of the apparent standard molar internal energy of a species.
    /// @param species The name of the species
    auto checkStandardInternalEnergy(std::string species) const -> bool;

    /// Return true if there is support for the calculation of the apparent standard molar enthalpy of a species.
    /// @param species The name of the species
    auto checkStandardEnthalpy(std::string species) const -> bool;

    /// Return true if there is support for the calculation of the standard molar entropies of a species.
    /// @param species The name of the species
    auto checkStandardEntropy(std::string species) const -> bool;

    /// Return true if there is support for the calculation of the standard molar volumes of a species.
    /// @param species The name of the species
    auto checkStandardVolume(std::string species) const -> bool;

    /// Return true if there is support for the calculation of the standard molar isobaric heat capacity of a species.
    /// @param species The name of the species
    auto checkStandardHeatCapacity(std::string species) const -> bool;

    /// Calculate the thermodynamic state of an aqueous species using the HKF model.
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    /// @param species The name of the species
    /// @see SpeciesThermoState
    auto speciesThermoStateHKF(double T, double P, std::string species) -> SpeciesThermoState;

    /// Calculate the thermodynamic state of water using the Haar-Gallagher-Kell (1984) equation of state.
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

} // namespace Reaktor
