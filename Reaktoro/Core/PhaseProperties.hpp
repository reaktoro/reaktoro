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
#include <Reaktoro/Common/ChemicalScalar.hpp>
#include <Reaktoro/Common/ChemicalVector.hpp>
#include <Reaktoro/Common/ThermoScalar.hpp>
#include <Reaktoro/Common/ThermoVector.hpp>

namespace Reaktoro {

// Forward declarations
class Phase;

/// Defines a class for querying thermodynamic properties of a phase.
class PhaseProperties
{
public:
    /// Construct a default PhaseProperties instance
    PhaseProperties();

    /// Return the temperature of the phase (in units of K).
    auto temperature() const -> ThermoScalar;

    /// Return the pressure of the phase (in units of Pa).
    auto pressure() const -> ThermoScalar;

    /// Return the amounts of the species of the phase (in units of mol).
    auto composition() const -> ChemicalVector;

    /// Return the molar fractions of the species.
    auto molarFractions() const -> ChemicalVector;

    /// Return the activity constants of the species.
    auto lnActivityConstants() const -> ChemicalVector;

    /// Return the activity coefficients of the species.
    auto lnActivityCoefficients() const -> ChemicalVector;

    /// Return the activities of the species.
    auto lnActivities() const -> ChemicalVector;

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

    /// Return the molar Gibbs energy of the phase (in units of J/mol).
    auto molarGibbsEnergy() const -> ChemicalScalar;

    /// Return the molar enthalpy of the phase (in units of J/mol).
    auto molarEnthalpy() const -> ChemicalScalar;

    /// Return the molar volume of the phase (in units of m3/mol).
    auto molarVolume() const -> ChemicalScalar;

    /// Return the molar entropy of the phase (in units of J/(mol*K)).
    auto molarEntropy() const -> ChemicalScalar;

    /// Return the molar internal energy of the phase (in units of J/mol).
    auto molarInternalEnergy() const -> ChemicalScalar;

    /// Return the molar Helmholtz energy of the phase (in units of J/mol).
    auto molarHelmholtzEnergy() const -> ChemicalScalar;

    /// Return the molar isobaric heat capacity of the phase (in units of J/(mol*K)).
    auto molarHeatCapacityConstP() const -> ChemicalScalar;

    /// Return the molar isochoric heat capacity of the phase (in units of J/(mol*K)).
    auto molarHeatCapacityConstV() const -> ChemicalScalar;

    /// Return the specific Gibbs energy of the phase (in units of J/kg).
    auto specificGibbsEnergy() const -> ChemicalScalar;

    /// Return the specific enthalpy of the phase (in units of J/kg).
    auto specificEnthalpy() const -> ChemicalScalar;

    /// Return the specific volume of the phase (in units of m3/kg).
    auto specificVolume() const -> ChemicalScalar;

    /// Return the specific entropy of the phase (in units of J/(kg*K)).
    auto specificEntropy() const -> ChemicalScalar;

    /// Return the specific internal energy of the phase (in units of J/kg).
    auto specificInternalEnergy() const -> ChemicalScalar;

    /// Return the specific Helmholtz energy of the phase (in units of J/kg).
    auto specificHelmholtzEnergy() const -> ChemicalScalar;

    /// Return the specific isobaric heat capacity of the phase (in units of J/(kg*K)).
    auto specificHeatCapacityConstP() const -> ChemicalScalar;

    /// Return the specific isochoric heat capacity of the phase (in units of J/(kg*K)).
    auto specificHeatCapacityConstV() const -> ChemicalScalar;

    /// Return the mass of the phase (in units of kg).
    auto mass() const -> ChemicalScalar;

    /// Return the number of moles in the phase (in units of mol).
    auto moles() const -> ChemicalScalar;

    /// Phase class is strongly coupled with PhaseProperties class
    friend class Phase;

private:
    /// The temperature of the phase (in units of K)
    ThermoScalar T;

    /// The pressure of the phase (in units of Pa)
    ThermoScalar P;

    /// The amounts of the species of the phase (in units of mol).
    ChemicalVector n;

    /// The molar fractions of the species of the phase (in units of mol).
    ChemicalVector x;

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

    /// The natural log of the activity constants of the species.
    ChemicalVector ln_activity_constants;

    /// The natural log of the activity coefficients of the species.
    ChemicalVector ln_activity_coefficients;

    /// The natural log of the activities of the species.
    ChemicalVector ln_activities;

    /// The molar Gibbs energy of the phase (in units of J/mol).
    ChemicalScalar molar_gibbs_energy;

    /// The molar enthalpy of the phase (in units of J/mol).
    ChemicalScalar molar_enthalpy;

    /// The molar volume of the phase (in units of m3/mol).
    ChemicalScalar molar_volume;

    /// The molar isobaric heat capacity of the phase (in units of J/(mol*K)).
    ChemicalScalar molar_heat_capacity_cp;

    /// The molar isochoric heat capacity of the phase (in units of J/(mol*K)).
    ChemicalScalar molar_heat_capacity_cv;

    /// The mass of the phase (in units of kg)
    ChemicalScalar total_mass;
};

} // namespace Reaktoro
