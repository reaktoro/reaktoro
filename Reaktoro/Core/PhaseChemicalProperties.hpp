//// Reaktoro is a unified framework for modeling chemically reactive systems.
////
//// Copyright (C) 2014-2015 Allan Leal
////
//// This program is free software: you can redistribute it and/or modify
//// it under the terms of the GNU General Public License as published by
//// the Free Software Foundation, either version 3 of the License, or
//// (at your option) any later version.
////
//// This program is distributed in the hope that it will be useful,
//// but WITHOUT ANY WARRANTY; without even the implied warranty of
//// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//// GNU General Public License for more details.
////
//// You should have received a copy of the GNU General Public License
//// along with this program. If not, see <http://www.gnu.org/licenses/>.
//
//#pragma once
//
//// C++ includes
//#include <memory>
//
//// Reaktoro includes
//#include <Reaktoro/Common/ChemicalVector.hpp>
//#include <Reaktoro/Common/ThermoVector.hpp>
//
//namespace Reaktoro {
//
//// Forward declarations
//class Phase;
//
///// Defines a class for querying thermodynamic and chemical properties of a phase.
//class PhaseChemicalProperties
//{
//public:
//    /// Construct a default PhaseChemicalProperties instance.
//    PhaseChemicalProperties();
//
//    /// Construct a PhaseChemicalProperties instance with given Phase instance.
//    explicit PhaseChemicalProperties(const Phase& phase);
//
//    /// Construct a copy of a PhaseChemicalProperties instance.
//    PhaseChemicalProperties(const PhaseChemicalProperties& other);
//
//    /// Destroy this instance.
//    virtual ~PhaseChemicalProperties();
//
//    /// Construct a copy of a PhaseChemicalProperties instance.
//    auto operator=(PhaseChemicalProperties other) -> PhaseChemicalProperties&;
//
//    /// Update the thermodynamic properties of the phase.
//    /// @param T The new temperature (in units of K)
//    /// @param P The new pressure (in units of Pa)
//    auto update(double T, double P) -> void;
//
//    /// Update the chemical properties of the phase.
//    /// @param T The new temperature (in units of K)
//    /// @param P The new pressure (in units of Pa)
//    /// @param n The new species amounts (in units of mol)
//    auto update(double T, double P, const Vector& n) -> void;
//
//    /// Return the temperature of the phase (in units of K).
//    auto temperature() const -> double;
//
//    /// Return the pressure of the phase (in units of Pa).
//    auto pressure() const -> double;
//
//    /// Return the amounts of the species of the phase (in units of mol).
//    auto composition() const -> const Vector&;
//
//    /// Return the molar fractions of the species.
//    auto molarFractions() const -> ChemicalVector;
//
//    /// Return the ln activity coefficients of the species.
//    auto lnActivityCoefficients() const -> ChemicalVector;
//
//    /// Return the ln activity constants of the species.
//    auto lnActivityConstants() const -> ThermoVector;
//
//    /// Return the ln activities of the species.
//    auto lnActivities() const -> ChemicalVector;
//
//    /// Return the chemical potentials of the species (in units of J/mol).
//    auto chemicalPotentials() const -> ChemicalVector;
//
//    /// Return the standard partial molar Gibbs energies of the species (in units of J/mol).
//    auto standardPartialMolarGibbsEnergies() const -> ThermoVector;
//
//    /// Return the standard partial molar enthalpies of the species (in units of J/mol).
//    auto standardPartialMolarEnthalpies() const -> ThermoVector;
//
//    /// Return the standard partial molar volumes of the species (in units of m3/mol).
//    auto standardPartialMolarVolumes() const -> ThermoVector;
//
//    /// Return the standard partial molar entropies of the species (in units of J/(mol*K)).
//    auto standardPartialMolarEntropies() const -> ThermoVector;
//
//    /// Return the standard partial molar internal energies of the species (in units of J/mol).
//    auto standardPartialMolarInternalEnergies() const -> ThermoVector;
//
//    /// Return the standard partial molar Helmholtz energies of the species (in units of J/mol).
//    auto standardPartialMolarHelmholtzEnergies() const -> ThermoVector;
//
//    /// Return the standard partial molar isobaric heat capacities of the species (in units of J/(mol*K)).
//    auto standardPartialMolarHeatCapacitiesConstP() const -> ThermoVector;
//
//    /// Return the standard partial molar isochoric heat capacities of the species (in units of J/(mol*K)).
//    auto standardPartialMolarHeatCapacitiesConstV() const -> ThermoVector;
//
//    /// Return the molar Gibbs energy of the phase (in units of J/mol).
//    auto molarGibbsEnergy() const -> ChemicalScalar;
//
//    /// Return the molar enthalpy of the phase (in units of J/mol).
//    auto molarEnthalpy() const -> ChemicalScalar;
//
//    /// Return the molar volume of the phase (in units of m3/mol).
//    auto molarVolume() const -> ChemicalScalar;
//
//    /// Return the molar entropy of the phase (in units of J/(mol*K)).
//    auto molarEntropy() const -> ChemicalScalar;
//
//    /// Return the molar internal energy of the phase (in units of J/mol).
//    auto molarInternalEnergy() const -> ChemicalScalar;
//
//    /// Return the molar Helmholtz energy of the phase (in units of J/mol).
//    auto molarHelmholtzEnergy() const -> ChemicalScalar;
//
//    /// Return the molar isobaric heat capacity of the phase (in units of J/(mol*K)).
//    auto molarHeatCapacityConstP() const -> ChemicalScalar;
//
//    /// Return the molar isochoric heat capacity of the phase (in units of J/(mol*K)).
//    auto molarHeatCapacityConstV() const -> ChemicalScalar;
//
//    /// Return the specific Gibbs energy of the phase (in units of J/kg).
//    auto specificGibbsEnergy() const -> ChemicalScalar;
//
//    /// Return the specific enthalpy of the phase (in units of J/kg).
//    auto specificEnthalpy() const -> ChemicalScalar;
//
//    /// Return the specific volume of the phase (in units of m3/kg).
//    auto specificVolume() const -> ChemicalScalar;
//
//    /// Return the specific entropy of the phase (in units of J/(kg*K)).
//    auto specificEntropy() const -> ChemicalScalar;
//
//    /// Return the specific internal energy of the phase (in units of J/kg).
//    auto specificInternalEnergy() const -> ChemicalScalar;
//
//    /// Return the specific Helmholtz energy of the phase (in units of J/kg).
//    auto specificHelmholtzEnergy() const -> ChemicalScalar;
//
//    /// Return the specific isobaric heat capacity of the phase (in units of J/(kg*K)).
//    auto specificHeatCapacityConstP() const -> ChemicalScalar;
//
//    /// Return the specific isochoric heat capacity of the phase (in units of J/(kg*K)).
//    auto specificHeatCapacityConstV() const -> ChemicalScalar;
//
//    /// Return the density of the phase (in units of kg/m3).
//    auto density() const -> ChemicalScalar;
//
//    /// Return the mass of the phase (in units of kg).
//    auto mass() const -> ChemicalScalar;
//
//    /// Return the amount of the phase (in units of mol).
//    auto amount() const -> ChemicalScalar;
//
//    /// Return the volume of the phase (in units of m3).
//    auto volume() const -> ChemicalScalar;
//
//private:
//    struct Impl;
//
//    std::unique_ptr<Impl> pimpl;
//};
//
//} // namespace Reaktoro
