// Reaktoro is a C++ library for computational reaction modelling.
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

// Reaktoro includes
#include <Reaktoro/Common/ThermoScalar.hpp>

namespace Reaktoro {

// Forward declarations
class Species;

/// Defines a class for querying standard thermodynamic properties of a species.
class SpeciesProperties
{
public:
    /// Construct a default SpeciesProperties instance
    SpeciesProperties();

    /// Return the temperature of the species (in units of K).
    auto temperature() const -> double;

    /// Return the pressure of the species (in units of Pa).
    auto pressure() const -> double;

    /// Return the standard partial molar Gibbs energy of the species (in units of J/mol).
    auto standardPartialMolarGibbsEnergy() const -> ThermoScalar;

    /// Return the standard partial molar enthalpy of the species (in units of J/mol).
    auto standardPartialMolarEnthalpy() const -> ThermoScalar;

    /// Return the standard partial molar volume of the species (in units of m3/mol).
    auto standardPartialMolarVolume() const -> ThermoScalar;

    /// Return the standard partial molar entropy of the species (in units of J/(mol*K)).
    auto standardPartialMolarEntropy() const -> ThermoScalar;

    /// Return the standard partial molar internal energy of the species (in units of J/mol).
    auto standardPartialMolarInternalEnergy() const -> ThermoScalar;

    /// Return the standard partial molar Helmholtz energy of the species (in units of J/mol).
    auto standardPartialMolarHelmholtzEnergy() const -> ThermoScalar;

    /// Return the standard partial molar isobaric heat capacity of the species (in units of J/(mol*K)).
    auto standardPartialMolarHeatCapacityConstP() const -> ThermoScalar;

    /// Return the standard partial molar isochoric heat capacity of the species (in units of J/(mol*K)).
    auto standardPartialMolarHeatCapacityConstV() const -> ThermoScalar;

    /// Species class is strongly coupled with SpeciesProperties class
    friend class Species;

private:
    /// The temperature of the phase (in units of K)
    ThermoScalar T;

    /// The pressure of the phase (in units of Pa)
    ThermoScalar P;

    /// The standard partial molar Gibbs energy of the species (in units of J/mol).
    ThermoScalar standard_partial_molar_gibbs_energy;

    /// The standard partial molar enthalpy of the species (in units of J/mol).
    ThermoScalar standard_partial_molar_enthalpy;

    /// The standard partial molar volume of the species (in units of m3/mol).
    ThermoScalar standard_partial_molar_volume;

    /// The standard partial molar isobaric heat capacity of the species (in units of J/(mol*K)).
    ThermoScalar standard_partial_molar_heat_capacity_cp;

    /// The standard partial molar isochoric heat capacity of the species (in units of J/(mol*K)).
    ThermoScalar standard_partial_molar_heat_capacity_cv;
};

} // namespace Reaktoro
