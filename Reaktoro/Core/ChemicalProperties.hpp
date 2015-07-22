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

/// Defines a class for querying thermodynamic and chemical properties of a set of phases.
class ChemicalProperties
{
public:
    /// Construct a default ChemicalProperties instance
    ChemicalProperties();

    /// Construct a ChemicalProperties instance with allocated memory
    /// @param nspecies The number of species in the system
    /// @param nphases The number of phases in the system
    ChemicalProperties(unsigned nspecies, unsigned nphases);

    /// Return the temperature of the phase (in units of K).
    auto temperature() const -> double;

    /// Return the pressure of the phase (in units of Pa).
    auto pressure() const -> double;

    /// Return the amounts of the species of the phase (in units of mol).
    auto composition() const -> Vector;

    /// Return the molar fractions of the species.
    auto molarFractions() const -> ChemicalVector;

    /// Return the activity coefficients of the species.
    auto lnActivityCoefficients() const -> ChemicalVector;

    /// Return the activities of the species.
    auto lnActivities() const -> ChemicalVector;

    /// Return the chemical potentials of the species (in units of J/mol).
    auto chemicalPotentials() const -> ChemicalVector;

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

    /// Return the molar Gibbs energies of the phases (in units of J/mol).
    auto phaseMolarGibbsEnergies() const -> ChemicalVector;

    /// Return the molar enthalpies of the phases (in units of J/mol).
    auto phaseMolarEnthalpies() const -> ChemicalVector;

    /// Return the molar volumes of the phases (in units of m3/mol).
    auto phaseMolarVolumes() const -> ChemicalVector;

    /// Return the molar entropies of the phases (in units of J/(mol*K)).
    auto phaseMolarEntropies() const -> ChemicalVector;

    /// Return the molar internal energies of the phases (in units of J/mol).
    auto phaseMolarInternalEnergies() const -> ChemicalVector;

    /// Return the molar Helmholtz energies of the phases (in units of J/mol).
    auto phaseMolarHelmholtzEnergies() const -> ChemicalVector;

    /// Return the molar isobaric heat capacities of the phases (in units of J/(mol*K)).
    auto phaseMolarHeatCapacitiesConstP() const -> ChemicalVector;

    /// Return the molar isochoric heat capacities of the phases (in units of J/(mol*K)).
    auto phaseMolarHeatCapacitiesConstV() const -> ChemicalVector;

    /// Return the specific Gibbs energies of the phases (in units of J/kg).
    auto phaseSpecificGibbsEnergies() const -> ChemicalVector;

    /// Return the specific enthalpies of the phases (in units of J/kg).
    auto phaseSpecificEnthalpies() const -> ChemicalVector;

    /// Return the specific volumes of the phases (in units of m3/kg).
    auto phaseSpecificVolumes() const -> ChemicalVector;

    /// Return the specific entropies of the phases (in units of J/(kg*K)).
    auto phaseSpecificEntropies() const -> ChemicalVector;

    /// Return the specific internal energies of the phases (in units of J/kg).
    auto phaseSpecificInternalEnergies() const -> ChemicalVector;

    /// Return the specific Helmholtz energies of the phases (in units of J/kg).
    auto phaseSpecificHelmholtzEnergies() const -> ChemicalVector;

    /// Return the specific isobaric heat capacities of the phases (in units of J/(kg*K)).
    auto phaseSpecificHeatCapacitiesConstP() const -> ChemicalVector;

    /// Return the specific isochoric heat capacities of the phases (in units of J/(kg*K)).
    auto phaseSpecificHeatCapacitiesConstV() const -> ChemicalVector;

    /// Return the masses of the phases (in units of kg).
    auto phaseMasses() const -> ChemicalVector;

    /// Return the molar amounts of the phases (in units of mol).
    auto phaseMoles() const -> ChemicalVector;

    /// Return the volumes of the phases (in units of m3).
    auto phaseVolumes() const -> ChemicalVector;

    // Class ChemicalSystem is strongly coupled with class ChemicalProperties
    friend class ChemicalSystem;

private:
    /// The temperature of the system (in units of K)
    ThermoScalar T;

    /// The pressure of the system (in units of Pa)
    ThermoScalar P;

    /// The amounts of the species in the system (in units of mol).
    ChemicalVector n;

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

    /// The molar fractions of the species of the phase (in units of mol).
    ChemicalVector molar_fractions;

    /// The natural log of the activity coefficients of the species.
    ChemicalVector ln_activity_coefficients;

    /// The natural log of the activities of the species.
    ChemicalVector ln_activities;

    /// The molar Gibbs energies of the phases (in units of J/mol).
    ChemicalVector phase_molar_gibbs_energies;

    /// The molar enthalpies of the phases (in units of J/mol).
    ChemicalVector phase_molar_enthalpies;

    /// The molar volumes of the phases (in units of m3/mol).
    ChemicalVector phase_molar_volumes;

    /// The molar isobaric heat capacities of the phases (in units of J/(mol*K)).
    ChemicalVector phase_molar_heat_capacities_cp;

    /// The molar isochoric heat capacities of the phases (in units of J/(mol*K)).
    ChemicalVector phase_molar_heat_capacities_cv;

    /// The mass of the phases (in units of mol)
    ChemicalVector phase_moles;

    /// The mass of the phases (in units of kg)
    ChemicalVector phase_masses;
};

/// Defines a class for querying thermodynamic and chemical properties of a phase.
class PhaseChemicalProperties
{
public:
    /// Construct a default PhaseChemicalProperties instance
    PhaseChemicalProperties();

    /// Construct a PhaseChemicalProperties instance with allocated memory
    /// @param nspecies The number of species in the phase
    explicit PhaseChemicalProperties(unsigned nspecies);

    /// Return the temperature of the phase (in units of K).
    auto temperature() const -> double;

    /// Return the pressure of the phase (in units of Pa).
    auto pressure() const -> double;

    /// Return the amounts of the species of the phase (in units of mol).
    auto composition() const -> Vector;

    /// Return the molar fractions of the species.
    auto molarFractions() const -> ChemicalVector;

    /// Return the activity coefficients of the species.
    auto lnActivityCoefficients() const -> ChemicalVector;

    /// Return the activities of the species.
    auto lnActivities() const -> ChemicalVector;

    /// Return the chemical potentials of the species (in units of J/mol).
    auto chemicalPotentials() const -> ChemicalVector;

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
    auto phaseMolarGibbsEnergy() const -> ChemicalScalar;

    /// Return the molar enthalpy of the phase (in units of J/mol).
    auto phaseMolarEnthalpy() const -> ChemicalScalar;

    /// Return the molar volume of the phase (in units of m3/mol).
    auto phaseMolarVolume() const -> ChemicalScalar;

    /// Return the molar entropy of the phase (in units of J/(mol*K)).
    auto phaseMolarEntropy() const -> ChemicalScalar;

    /// Return the molar internal energy of the phase (in units of J/mol).
    auto phaseMolarInternalEnergy() const -> ChemicalScalar;

    /// Return the molar Helmholtz energy of the phase (in units of J/mol).
    auto phaseMolarHelmholtzEnergy() const -> ChemicalScalar;

    /// Return the molar isobaric heat capacity of the phase (in units of J/(mol*K)).
    auto phaseMolarHeatCapacityConstP() const -> ChemicalScalar;

    /// Return the molar isochoric heat capacity of the phase (in units of J/(mol*K)).
    auto phaseMolarHeatCapacityConstV() const -> ChemicalScalar;

    /// Return the specific Gibbs energy of the phase (in units of J/kg).
    auto phaseSpecificGibbsEnergy() const -> ChemicalScalar;

    /// Return the specific enthalpy of the phase (in units of J/kg).
    auto phaseSpecificEnthalpy() const -> ChemicalScalar;

    /// Return the specific volume of the phase (in units of m3/kg).
    auto phaseSpecificVolume() const -> ChemicalScalar;

    /// Return the specific entropy of the phase (in units of J/(kg*K)).
    auto phaseSpecificEntropy() const -> ChemicalScalar;

    /// Return the specific internal energy of the phase (in units of J/kg).
    auto phaseSpecificInternalEnergy() const -> ChemicalScalar;

    /// Return the specific Helmholtz energy of the phase (in units of J/kg).
    auto phaseSpecificHelmholtzEnergy() const -> ChemicalScalar;

    /// Return the specific isobaric heat capacity of the phase (in units of J/(kg*K)).
    auto phaseSpecificHeatCapacityConstP() const -> ChemicalScalar;

    /// Return the specific isochoric heat capacity of the phase (in units of J/(kg*K)).
    auto phaseSpecificHeatCapacityConstV() const -> ChemicalScalar;

    /// Return the mass of the phase (in units of kg).
    auto phaseMass() const -> ChemicalScalar;

    /// Return the molar amounts of the phase (in units of mol).
    auto phaseMoles() const -> ChemicalScalar;

    /// Return the volumes of the phase (in units of m3).
    auto phaseVolume() const -> ChemicalScalar;

    // Class Phase is strongly coupled with class PhaseChemicalProperties
    friend class Phase;

private:
    /// The temperature of the phase (in units of K)
    ThermoScalar T;

    /// The pressure of the phase (in units of Pa)
    ThermoScalar P;

    /// The amounts of the species of the phase (in units of mol).
    ChemicalVector n;

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

    /// The molar fractions of the species of the phase (in units of mol).
    ChemicalVector molar_fractions;

    /// The natural log of the activity coefficients of the species.
    ChemicalVector ln_activity_coefficients;

    /// The natural log of the activities of the species.
    ChemicalVector ln_activities;

    /// The molar Gibbs energy of the phase (in units of J/mol).
    ChemicalScalar phase_molar_gibbs_energy;

    /// The molar enthalpy of the phase (in units of J/mol).
    ChemicalScalar phase_molar_enthalpy;

    /// The molar volumes of the phase (in units of m3/mol).
    ChemicalScalar phase_molar_volume;

    /// The molar isobaric heat capacity of the phase (in units of J/(mol*K)).
    ChemicalScalar phase_molar_heat_capacity_cp;

    /// The molar isochoric heat capacity of the phase (in units of J/(mol*K)).
    ChemicalScalar phase_molar_heat_capacity_cv;

    /// The molar amount of the phase (in units of mol)
    ChemicalScalar phase_moles;

    /// The mass of the phase (in units of kg)
    ChemicalScalar phase_mass;
};

} // namespace Reaktoro
