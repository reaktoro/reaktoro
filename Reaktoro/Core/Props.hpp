// Reaktoro is a unified framework for modeling chemically reactive phases.
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

// Reaktoro includes
#include <Reaktoro/Common/ChemicalScalar.hpp>
#include <Reaktoro/Common/ChemicalVector.hpp>
#include <Reaktoro/Common/ThermoVector.hpp>
#include <Reaktoro/Core/ActivityModel.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Core/StandardThermoModel.hpp>

namespace Reaktoro {

/// The chemical properties of a phase and its chemical species.
class PhaseChemicalProps
{
public:
    /// Construct a default PhaseChemicalProps instance.
    PhaseChemicalProps();

    /// Return the temperature of the phase (in units of K).
    auto temperature() const -> Temperature;

    /// Return the pressure of the phase (in units of Pa).
    auto pressure() const -> Pressure;

    /// Return the amounts of the species in the phase (in units of mol).
    auto speciesAmounts() const -> Composition;

    /// Return the mole fractions of the species in the phase.
    auto moleFractions() const -> ChemicalVector;

    /// Return the ln activity coefficients of the species in the phase.
    auto lnActivityCoefficients() const -> ChemicalVectorConstRef;

    /// Return the ln activities of the species in the phase.
    auto lnActivities() const -> ChemicalVectorConstRef;

    /// Return the partial molar volume of the species in the phase (in units of m3/mol).
    auto partialMolarVolumes() const -> ChemicalVectorConstRef;

    /// Return the chemical potentials of the species (in units of J/mol).
    auto chemicalPotentials() const -> ChemicalVector;

    /// Return the standard partial molar Gibbs energies of the species (in units of J/mol).
    auto standardPartialMolarGibbsEnergies() const -> ThermoVectorConstRef;

    /// Return the standard partial molar enthalpies of the species (in units of J/mol).
    auto standardPartialMolarEnthalpies() const -> ThermoVectorConstRef;

    /// Return the standard partial molar volumes of the species (in units of m3/mol).
    auto standardPartialMolarVolumes() const -> ThermoVectorConstRef;

    /// Return the standard partial molar entropies of the species (in units of J/(mol*K)).
    auto standardPartialMolarEntropies() const -> ThermoVector;

    /// Return the standard partial molar internal energies of the species (in units of J/mol).
    auto standardPartialMolarInternalEnergies() const -> ThermoVector;

    /// Return the standard partial molar Helmholtz energies of the species (in units of J/mol).
    auto standardPartialMolarHelmholtzEnergies() const -> ThermoVector;

    /// Return the standard partial molar isobaric heat capacities of the species (in units of J/(mol*K)).
    auto standardPartialMolarHeatCapacitiesConstP() const -> ThermoVectorConstRef;

    /// Return the standard partial molar isochoric heat capacities of the species (in units of J/(mol*K)).
    auto standardPartialMolarHeatCapacitiesConstV() const -> ThermoVectorConstRef;

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
    auto specificVolumes() const -> ChemicalScalar;

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

    /// Return the density of the phase (in units of kg/m3).
    auto density() const -> ChemicalScalar;

    /// Return the mass of the phase (in units of kg).
    auto mass() const -> ChemicalScalar;

    /// Return the total amount of species in the phase (in units of mol).
    auto amount() const -> ChemicalScalar;

    /// Return the volume of the phase (in units of m3).
    auto volume() const -> ChemicalScalar;

private:
    /// The number of species in the phase
    Index num_species;

    /// The temperature of the phase (in units of K)
    Temperature T;

    /// The pressure of the phase (in units of Pa)
    Pressure P;

    /// The amounts of the species in the phase (in units of mol).
    Vector n;

    /// The mole fractions of the species in the phase (in units of mol/mol).
    ChemicalVector x;

    /// The standard thermodynamic properties of each species in the phase.
    std::vector<StandardThermoProps> species_standard_thermo_props;

    /// The activity and excess thermodynamic properties of the phase and its species.
    ActivityProps activity_props;

    // IMPORTANT: Needed for performance reasons.
    friend Phase;
};

} // namespace Reaktoro
