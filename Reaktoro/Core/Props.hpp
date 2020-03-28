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
#include <Reaktoro/Common/Real.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Core/ActivityModel.hpp>
#include <Reaktoro/Core/StandardThermoModel.hpp>
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

/// The chemical properties of a phase and its chemical species.
class PhaseChemicalProps
{
public:
    /// Construct a default PhaseChemicalProps instance.
    PhaseChemicalProps();

    /// Return the temperature of the phase (in units of K).
    auto temperature() const -> real;

    /// Return the pressure of the phase (in units of Pa).
    auto pressure() const -> real;

    /// Return the amounts of the species in the phase (in units of mol).
    auto speciesAmounts() const -> VectorXrConstRef;

    /// Return the mole fractions of the species in the phase.
    auto moleFractions() const -> VectorXd;

    /// Return the ln activity coefficients of the species in the phase.
    auto lnActivityCoefficients() const -> VectorXdConstRef;

    /// Return the ln activities of the species in the phase.
    auto lnActivities() const -> VectorXdConstRef;

    /// Return the partial molar volume of the species in the phase (in units of m3/mol).
    auto partialMolarVolumes() const -> VectorXdConstRef;

    /// Return the chemical potentials of the species (in units of J/mol).
    auto chemicalPotentials() const -> VectorXd;

    /// Return the standard partial molar Gibbs energies of the species (in units of J/mol).
    auto standardPartialMolarGibbsEnergies() const -> VectorXrConstRef;

    /// Return the standard partial molar enthalpies of the species (in units of J/mol).
    auto standardPartialMolarEnthalpies() const -> VectorXrConstRef;

    /// Return the standard partial molar volumes of the species (in units of m3/mol).
    auto standardPartialMolarVolumes() const -> VectorXrConstRef;

    /// Return the standard partial molar entropies of the species (in units of J/(mol*K)).
    auto standardPartialMolarEntropies() const -> VectorXr;

    /// Return the standard partial molar internal energies of the species (in units of J/mol).
    auto standardPartialMolarInternalEnergies() const -> VectorXr;

    /// Return the standard partial molar Helmholtz energies of the species (in units of J/mol).
    auto standardPartialMolarHelmholtzEnergies() const -> VectorXr;

    /// Return the standard partial molar isobaric heat capacities of the species (in units of J/(mol*K)).
    auto standardPartialMolarHeatCapacitiesConstP() const -> VectorXrConstRef;

    /// Return the standard partial molar isochoric heat capacities of the species (in units of J/(mol*K)).
    auto standardPartialMolarHeatCapacitiesConstV() const -> VectorXrConstRef;

    /// Return the molar Gibbs energy of the phase (in units of J/mol).
    auto molarGibbsEnergy() const -> real;

    /// Return the molar enthalpy of the phase (in units of J/mol).
    auto molarEnthalpy() const -> real;

    /// Return the molar volume of the phase (in units of m3/mol).
    auto molarVolume() const -> real;

    /// Return the molar entropy of the phase (in units of J/(mol*K)).
    auto molarEntropy() const -> real;

    /// Return the molar internal energy of the phase (in units of J/mol).
    auto molarInternalEnergy() const -> real;

    /// Return the molar Helmholtz energy of the phase (in units of J/mol).
    auto molarHelmholtzEnergy() const -> real;

    /// Return the molar isobaric heat capacity of the phase (in units of J/(mol*K)).
    auto molarHeatCapacityConstP() const -> real;

    /// Return the molar isochoric heat capacity of the phase (in units of J/(mol*K)).
    auto molarHeatCapacityConstV() const -> real;

    /// Return the specific Gibbs energy of the phase (in units of J/kg).
    auto specificGibbsEnergy() const -> real;

    /// Return the specific enthalpy of the phase (in units of J/kg).
    auto specificEnthalpy() const -> real;

    /// Return the specific volume of the phase (in units of m3/kg).
    auto specificVolumes() const -> real;

    /// Return the specific entropy of the phase (in units of J/(kg*K)).
    auto specificEntropy() const -> real;

    /// Return the specific internal energy of the phase (in units of J/kg).
    auto specificInternalEnergy() const -> real;

    /// Return the specific Helmholtz energy of the phase (in units of J/kg).
    auto specificHelmholtzEnergy() const -> real;

    /// Return the specific isobaric heat capacity of the phase (in units of J/(kg*K)).
    auto specificHeatCapacityConstP() const -> real;

    /// Return the specific isochoric heat capacity of the phase (in units of J/(kg*K)).
    auto specificHeatCapacityConstV() const -> real;

    /// Return the density of the phase (in units of kg/m3).
    auto density() const -> real;

    /// Return the mass of the phase (in units of kg).
    auto mass() const -> real;

    /// Return the total amount of species in the phase (in units of mol).
    auto amount() const -> real;

    /// Return the volume of the phase (in units of m3).
    auto volume() const -> real;

private:
    /// The number of species in the phase
    Index num_species;

    /// The temperature of the phase (in units of K)
    real T;

    /// The pressure of the phase (in units of Pa)
    real P;

    /// The amounts of the species in the phase (in units of mol).
    VectorXr n;

    /// The mole fractions of the species in the phase (in units of mol/mol).
    VectorXd x;

    /// The standard thermodynamic properties of each species in the phase.
    std::vector<StandardThermoProps> species_standard_thermo_props;

    /// The activity and excess thermodynamic properties of the phase and its species.
    ActivityProps activity_props;

    // IMPORTANT: Needed for performance reasons.
    friend Phase;
};

} // namespace Reaktoro
