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
#include <Reaktoro/Common/ThermoVector.hpp>

namespace Reaktoro {

/// Defines a class for querying standard thermodynamic properties of a set of species.
class ThermoProperties
{
public:
    /// Construct a default ThermoProperties instance
    ThermoProperties();

    /// Construct a ThermoProperties instance with allocated memory
    explicit ThermoProperties(unsigned nspecies);

    /// Return the temperature of the phase (in units of K).
    auto temperature() const -> Temperature;

    /// Return the pressure of the phase (in units of Pa).
    auto pressure() const -> Pressure;

    /// Return the standard partial molar Gibbs energies of the species (in units of J/mol).
    auto standardPartialMolarGibbsEnergies() const -> ThermoVector;

    /// Return the standard partial molar enthalpies of the species (in units of J/mol).
    auto standardPartialMolarEnthalpies() const -> ThermoVector;

    /// Return the standard partial molar volumes of the species (in units of m3/mol).
    auto standardPartialMolarVolumes() const -> ThermoVector;

    /// Return the standard partial molar entropies of the species (in units of J/(mol*K)).
    auto standardPartialMolarEntropies() const -> ThermoVector;

    /// Return the standard partial molar internal energies of the species (in units of J/mol).
    auto standardPartialMolarInternalEnergies() const -> ThermoVector;

    /// Return the standard partial molar Helmholtz energies of the species (in units of J/mol).
    auto standardPartialMolarHelmholtzEnergies() const -> ThermoVector;

    /// Return the standard partial molar isobaric heat capacities of the species (in units of J/(mol*K)).
    auto standardPartialMolarHeatCapacitiesConstP() const -> ThermoVector;

    /// Return the standard partial molar isochoric heat capacities of the species (in units of J/(mol*K)).
    auto standardPartialMolarHeatCapacitiesConstV() const -> ThermoVector;

    // Both classes Phase and ChemicalSystem are strongly coupled with class ThermoProperties
    friend class Phase;
    friend class ChemicalSystem;

private:
    /// The temperature of the phase (in units of K)
    Temperature T;

    /// The pressure of the phase (in units of Pa)
    Pressure P;

    /// The standard partial molar Gibbs energies of the species (in units of J/mol).
    ThermoVector standard_partial_molar_gibbs_energies;

    /// The standard partial molar enthalpies of the species (in units of J/mol).
    ThermoVector standard_partial_molar_enthalpies;

    /// The standard partial molar volumes of the species (in units of m3/mol).
    ThermoVector standard_partial_molar_volumes;

    /// The standard partial molar isobaric heat capacities of the species (in units of J/(mol*K)).
    ThermoVector standard_partial_molar_heat_capacities_cp;

    /// The standard partial molar isochoric heat capacities of the species (in units of J/(mol*K)).
    ThermoVector standard_partial_molar_heat_capacities_cv;
};

/// Defines a class for querying standard thermodynamic properties of a species.
class SpeciesThermoProperties
{
public:
    /// Construct a default SpeciesThermoProperties instance
    SpeciesThermoProperties();

    /// Return the temperature of the species (in units of K).
    auto temperature() const -> Temperature;

    /// Return the pressure of the species (in units of Pa).
    auto pressure() const -> Pressure;

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

    // Class Species is strongly coupled with class SpeciesThermoProperties
    friend class Species;

private:
    /// The temperature of the phase (in units of K)
    Temperature T;

    /// The pressure of the phase (in units of Pa)
    Pressure P;

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
