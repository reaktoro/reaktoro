// Reaktoro is a unified framework for modeling chemically reactive phases.
//
// Copyright © 2014-2021 Allan Leal
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
#include <Reaktoro/Common/ArrayStream.hpp>
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/ChemicalPropsPhase.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalState;

/// The class that computes chemical properties of a chemical system.
class ChemicalProps
{
public:
    /// Construct an uninitialized ChemicalProps object with given chemical system.
    explicit ChemicalProps(const ChemicalSystem& system);

    /// Construct a ChemicalProps object with given chemical state of the system.
    explicit ChemicalProps(const ChemicalState& state);

    /// Update the chemical properties of the system.
    /// @param state The chemical state of the system
    auto update(const ChemicalState& state) -> void;

    /// Update the chemical properties of the system.
    /// @param T The temperature condition (in K)
    /// @param P The pressure condition (in Pa)
    /// @param n The amounts of the species in the system (in mol)
    auto update(const real& T, const real& P, ArrayXrConstRef n) -> void;

    /// Update the chemical properties of the system with serialized data.
    /// @param u The chemical properties of the system serialized in an array of real numbers.
    auto update(ArrayXrConstRef u) -> void;

    /// Update the chemical properties of the system with serialized data.
    /// @param u The chemical properties of the system serialized in an array of double numbers.
    auto update(ArrayXdConstRef u) -> void;

    /// Update the chemical properties of the system using ideal activity models.
    /// @param state The chemical state of the system
    auto updateIdeal(const ChemicalState& state) -> void;

    /// Update the chemical properties of the system using ideal activity models.
    /// @param T The temperature condition (in K)
    /// @param P The pressure condition (in Pa)
    /// @param n The amounts of the species in the system (in mol)
    auto updateIdeal(const real& T, const real& P, ArrayXrConstRef n) -> void;

    /// Serialize the chemical properties into the array stream @p stream.
    /// @param stream The array stream used to serialize the chemical properties.
    auto serialize(ArrayStream<real>& stream) const -> void;

    /// Serialize the chemical properties into the array stream @p stream.
    /// @param stream The array stream used to serialize the chemical properties.
    auto serialize(ArrayStream<double>& stream) const -> void;

    /// Update the chemical properties of the system using the array stream @p stream.
    /// @param stream The array stream containing the serialized chemical properties.
    auto deserialize(const ArrayStream<real>& stream) -> void;

    /// Update the chemical properties of the system using the array stream @p stream.
    /// @param stream The array stream containing the serialized chemical properties.
    auto deserialize(const ArrayStream<double>& stream) -> void;

    /// Return the chemical system associated with these chemical properties.
    auto system() const -> const ChemicalSystem&;

    /// Return the chemical properties of a phase with given index.
    /// @param phase The name or index of the phase in the system.
    auto phaseProps(StringOrIndex phase) -> ChemicalPropsPhaseRef;

    /// Return the chemical properties of a phase with given index.
    /// @param phase The name or index of the phase in the system.
    auto phaseProps(StringOrIndex phase) const -> ChemicalPropsPhaseConstRef;

    /// Return the extra data produced during the evaluation of activity models.
    auto extra() const -> const Map<String, Any>&;

    /// Return the temperature of the system (in K).
    auto temperature() const -> real;

    /// Return the pressure of the system (in Pa).
    auto pressure() const -> real;

    /// Return the amount of electric charge in the system (in mol).
    auto charge() const -> real;

    /// Return the amount of an element in the system (in mol).
    /// @param element The symbol or index of the element in the system.
    auto elementAmount(StringOrIndex element) const -> real;

    /// Return the amount of an element in the system (in mol).
    /// @param element The symbol or index of the element in the system.
    /// @param phase The name or index of the phase in the system.
    auto elementAmountInPhase(StringOrIndex element, StringOrIndex phase) const -> real;

    /// Return the amount of an element among a group of species in the system (in mol).
    /// @param element The symbol or index of the element in the system.
    /// @param indices The indices of the species in the system.
    auto elementAmountAmongSpecies(StringOrIndex element, ArrayXlConstRef indices) const -> real;

    /// Return the mass of an element in the system (in kg).
    /// @param element The symbol or index of the element in the system.
    auto elementMass(StringOrIndex element) const -> real;

    /// Return the mass of an element in the system (in kg).
    /// @param element The symbol or index of the element in the system.
    /// @param phase The name or index of the phase in the system.
    auto elementMassInPhase(StringOrIndex element, StringOrIndex phase) const -> real;

    /// Return the mass of an element among a group of species in the system (in kg).
    /// @param element The symbol or index of the element in the system.
    /// @param indices The indices of the species in the system.
    auto elementMassAmongSpecies(StringOrIndex element, ArrayXlConstRef indices) const -> real;

    /// Return the amounts of the elements in the system (in mol).
    auto elementAmounts() const -> ArrayXr;

    /// Return the amounts of the elements in a phase of the system (in mol).
    /// @param phase The name or index of the phase in the system.
    auto elementAmountsInPhase(StringOrIndex phase) const -> ArrayXr;

    /// Return the amounts of the elements among a group of species in the system (in mol).
    /// @param indices The indices of the species in the system.
    auto elementAmountsAmongSpecies(ArrayXlConstRef indices) const -> ArrayXr;

    /// Return the amount of a species in the system (in mol).
    /// @param species The name or index of the species in the system.
    auto speciesAmount(StringOrIndex species) const -> real;

    /// Return the mass of a species in the system.
    /// @param species The name or index of the species in the system.
    auto speciesMass(StringOrIndex species) const -> real;

    /// Return the mole fraction of a species in the system.
    /// @param species The name or index of the species in the system.
    auto speciesMoleFraction(StringOrIndex species) const -> real;

    /// Return the concentration (activity divided by activity coefficient) of a species in the system.
    /// @param species The name or index of the species in the system.
    auto speciesConcentration(StringOrIndex species) const -> real;

    /// Return the lg concentration (activity divided by activity coefficient) of a species in the system.
    /// @param species The name or index of the species in the system.
    auto speciesConcentrationLg(StringOrIndex species) const -> real;

    /// Return the ln concentration (activity divided by activity coefficient) of a species in the system.
    /// @param species The name or index of the species in the system.
    auto speciesConcentrationLn(StringOrIndex species) const -> real;

    /// Return the activity coefficient of a species in the system.
    /// @param species The name or index of the species in the system.
    auto speciesActivityCoefficient(StringOrIndex species) const -> real;

    /// Return the lg activity coefficient of a species in the system.
    /// @param species The name or index of the species in the system.
    auto speciesActivityCoefficientLg(StringOrIndex species) const -> real;

    /// Return the ln activity coefficient of a species in the system.
    /// @param species The name or index of the species in the system.
    auto speciesActivityCoefficientLn(StringOrIndex species) const -> real;

    /// Return the activity of a species in the system.
    /// @param species The name or index of the species in the system.
    auto speciesActivity(StringOrIndex species) const -> real;

    /// Return the lg activity of a species in the system.
    /// @param species The name or index of the species in the system.
    auto speciesActivityLg(StringOrIndex species) const -> real;

    /// Return the ln activity of a species in the system.
    /// @param species The name or index of the species in the system.
    auto speciesActivityLn(StringOrIndex species) const -> real;

    /// Return the chemical potential of a species in the system.
    /// @param species The name or index of the species in the system.
    auto speciesChemicalPotential(StringOrIndex species) const -> real;

    /// Return the standard partial molar volume of a species in the system (in m³/mol).
    /// @param species The name or index of the species in the system.
    auto speciesStandardVolume(StringOrIndex species) const -> real;

    /// Return the temperature derivative of the standard partial molar volume of a species in the system (in m³/(mol·K)).
    /// @param species The name or index of the species in the system.
    auto speciesStandardVolumeT(StringOrIndex species) const -> real;

    /// Return the pressure derivative of the standard partial molar volume of a species in the system (in m³/(mol·Pa)).
    /// @param species The name or index of the species in the system.
    auto speciesStandardVolumeP(StringOrIndex species) const -> real;

    /// Return the standard partial molar Gibbs energy of formation of a species in the system (in J/mol).
    /// @param species The name or index of the species in the system.
    auto speciesStandardGibbsEnergy(StringOrIndex species) const -> real;

    /// Return the standard partial molar enthalpy of formation of a species in the system (in J/mol).
    /// @param species The name or index of the species in the system.
    auto speciesStandardEnthalpy(StringOrIndex species) const -> real;

    /// Return the standard partial molar entropy of formation of the species a the system (in J/(mol·K)).
    /// @param species The name or index of the species in the system.
    auto speciesStandardEntropy(StringOrIndex species) const -> real;

    /// Return the standard partial molar internal energy of formation of a species in the system (in J/mol).
    /// @param species The name or index of the species in the system.
    auto speciesStandardInternalEnergy(StringOrIndex species) const -> real;

    /// Return the standard partial molar Helmholtz energy of formation of a species in the system (in J/mol).
    /// @param species The name or index of the species in the system.
    auto speciesStandardHelmholtzEnergy(StringOrIndex species) const -> real;

    /// Return the standard partial molar isobaric heat capacity of the species a the system (in J/(mol·K)).
    /// @param species The name or index of the species in the system.
    auto speciesStandardHeatCapacityConstP(StringOrIndex species) const -> real;

    /// Return the standard partial molar isochoric heat capacity of the species a the system (in J/(mol·K)).
    /// @param species The name or index of the species in the system.
    auto speciesStandardHeatCapacityConstV(StringOrIndex species) const -> real;

    /// Return the amounts of the species in the system (in mol).
    auto speciesAmounts() const -> ArrayXrConstRef;

    /// Return the masses of the species in the system (in kg).
    auto speciesMasses() const -> ArrayXr;

    /// Return the mole fractions of the species in the system.
    auto speciesMoleFractions() const -> ArrayXrConstRef;

    /// Return the ln concentrations (activity divided by activity coefficient) of the species in the system.
    auto speciesConcentrationsLn() const -> ArrayXr;

    /// Return the ln activity coefficients of the species in the system.
    auto speciesActivityCoefficientsLn() const -> ArrayXrConstRef;

    /// Return the ln activities of the species in the system.
    auto speciesActivitiesLn() const -> ArrayXrConstRef;

    /// Return the chemical potentials of the species in the system (in J/mol).
    auto speciesChemicalPotentials() const -> ArrayXrConstRef;

    /// Return the standard partial molar volumes of the species in the system (in m³/mol).
    auto speciesStandardVolumes() const -> ArrayXrConstRef;

    /// Return the temperature derivative of the standard molar volumes of the species in the system (in m³/(mol·K)).
    auto speciesStandardVolumesT() const -> ArrayXrConstRef;

    /// Return the pressure derivative of the standard molar volumes of the species in the system (in m³/(mol·Pa)).
    auto speciesStandardVolumesP() const -> ArrayXrConstRef;

    /// Return the standard partial molar Gibbs energies of formation of the species in the system (in J/mol).
    auto speciesStandardGibbsEnergies() const -> ArrayXrConstRef;

    /// Return the standard partial molar enthalpies of formation of the species in the system (in J/mol).
    auto speciesStandardEnthalpies() const -> ArrayXrConstRef;

    /// Return the standard partial molar entropies of formation of the species in the system (in J/(mol·K)).
    auto speciesStandardEntropies() const -> ArrayXr;

    /// Return the standard partial molar internal energies of formation of the species in the system (in J/mol).
    auto speciesStandardInternalEnergies() const -> ArrayXr;

    /// Return the standard partial molar Helmholtz energies of formation of the species in the system (in J/mol).
    auto speciesStandardHelmholtzEnergies() const -> ArrayXr;

    /// Return the standard partial molar isobaric heat capacities of the species in the system (in J/(mol·K)).
    auto speciesStandardHeatCapacitiesConstP() const -> ArrayXrConstRef;

    /// Return the standard partial molar isochoric heat capacities of the species in the system (in J/(mol·K)).
    auto speciesStandardHeatCapacitiesConstV() const -> ArrayXr;

    /// Return the molar volume of the system (in m³/mol).
    auto molarVolume() const -> real;

    /// Return the temperature derivative of the molar volume of the system (in m³/(mol·K)).
    auto molarVolumeT() const -> real;

    /// Return the pressure derivative of the molar volume of the system (in m³/(mol·Pa)).
    auto molarVolumeP() const -> real;

    /// Return the molar Gibbs energy of formation of the system (in J/mol).
    auto molarGibbsEnergy() const -> real;

    /// Return the molar enthalpy of formation of the system (in J/mol).
    auto molarEnthalpy() const -> real;

    /// Return the molar entropy of formation of the system (in J/(mol·K)).
    auto molarEntropy() const -> real;

    /// Return the molar internal energy of formation of the system (in J/mol).
    auto molarInternalEnergy() const -> real;

    /// Return the molar Helmholtz energy of formation of the system (in J/mol).
    auto molarHelmholtzEnergy() const -> real;

    /// Return the molar isobaric heat capacity of the system (in J/(mol·K)).
    auto molarHeatCapacityConstP() const -> real;

    /// Return the molar isochoric heat capacity of the system (in J/(mol·K)).
    auto molarHeatCapacityConstV() const -> real;

    /// Return the specific volume of the system (in m³/kg).
    auto specificVolume() const -> real;

    /// Return the temperature derivative of the specific volume of the system (in m³/(kg·K)).
    auto specificVolumeT() const -> real;

    /// Return the pressure derivative of the specific volume of the system (in m³/(kg·Pa)).
    auto specificVolumeP() const -> real;

    /// Return the specific Gibbs energy of formation of the system (in J/kg).
    auto specificGibbsEnergy() const -> real;

    /// Return the specific enthalpy of formation of the system (in J/kg).
    auto specificEnthalpy() const -> real;

    /// Return the specific entropy of formation of the system (in J/(kg·K)).
    auto specificEntropy() const -> real;

    /// Return the specific internal energy of formation of the system (in J/kg).
    auto specificInternalEnergy() const -> real;

    /// Return the specific Helmholtz energy of formation of the system (in J/kg).
    auto specificHelmholtzEnergy() const -> real;

    /// Return the specific isobaric heat capacity of the system (in J/(kg·K)).
    auto specificHeatCapacityConstP() const -> real;

    /// Return the specific isochoric heat capacity of the system (in J/(kg·K)).
    auto specificHeatCapacityConstV() const -> real;

    /// Return the density of the system (in kg/m³).
    auto density() const -> real;

    /// Return the sum of species amounts in the system (in mol).
    auto amount() const -> real;

    /// Return the sum of species masses in the system (in kg).
    auto mass() const -> real;

    /// Return the volume of the system (in m³).
    auto volume() const -> real;

    /// Return the temperature derivative of the volume of the system (in m³/K).
    auto volumeT() const -> real;

    /// Return the pressure derivative of the volume of the system (in m³/Pa).
    auto volumeP() const -> real;

    /// Return the Gibbs energy of formation of the system (in J).
    auto gibbsEnergy() const -> real;

    /// Return the enthalpy of formation of the system (in J).
    auto enthalpy() const -> real;

    /// Return the entropy of formation of the system (in J/K).
    auto entropy() const -> real;

    /// Return the internal energy of formation of the system (in J).
    auto internalEnergy() const -> real;

    /// Return the Helmholtz energy of formation of the system (in J).
    auto helmholtzEnergy() const -> real;

    /// Return the isobaric heat capacity of the system (in J/K).
    auto heatCapacityConstP() const -> real;

    /// Return the isochoric heat capacity of the system (in J/K).
    auto heatCapacityConstV() const -> real;

    /// Output the chemical properties of the system to a stream.
    auto output(std::ostream& out) const -> void;

    /// Output the chemical properties of the system to a file.
    auto output(const String& filename) const -> void;

    /// Return the chemical properties in this object serialized in an array of real numbers.
    operator VectorXr() const;

    /// Return the chemical properties in this object serialized in an array of double numbers.
    operator VectorXd() const;

private:
    /// The chemical system associated with these chemical properties.
    ChemicalSystem msystem;

    /// The temperature of the system (in K).
    real T = 0.0;

    /// The pressure of the system (in Pa).
    real P = 0.0;

    /// The temperatures of each phase (in K).
    ArrayXr Ts;

    /// The pressures of each phase (in Pa).
    ArrayXr Ps;

    /// The amounts of each species in the system (in mol).
    ArrayXr n;

    /// The sum of species amounts in each phase of the system (in mol).
    ArrayXr nsum;

    /// The sum of species masses in each phase of the system (in kg).
    ArrayXr msum;

    /// The mole fractions of the species in the system (in mol/mol).
    ArrayXr x;

    /// The standard molar Gibbs energies of formation of the species in the system (in J/mol).
    ArrayXr G0;

    /// The standard molar enthalpies of formation of the species in the system (in J/mol).
    ArrayXr H0;

    /// The standard molar volumes of the species in the system (in m³/mol).
    ArrayXr V0;

    /// The temperature derivative of the standard molar volumes of the species in the system (in m³/(mol·K)).
    ArrayXr VT0;

    /// The pressure derivative of the standard molar volumes of the species in the system (in m³/(mol·Pa)).
    ArrayXr VP0;

    /// The standard molar isobaric heat capacities of the species in the system (in J/(mol·K)).
    ArrayXr Cp0;

    /// The corrective molar volume of each phase in the system (in m³/mol).
    ArrayXr Vx;

    /// The temperature derivative at constant pressure of the corrective molar volume of each phase in the system (in m³/(mol·K)).
    ArrayXr VxT;

    /// The pressure derivative at constant temperature of the corrective molar volume of each phase in the system (in m³/(mol·Pa)).
    ArrayXr VxP;

    /// The corrective molar Gibbs energy of each phase in the system (in J/mol).
    ArrayXr Gx;

    /// The corrective molar enthalpy of each phase in the system (in J/mol).
    ArrayXr Hx;

    /// The corrective molar isobaric heat capacity of each phase in the system (in J/(mol·K)).
    ArrayXr Cpx;

    /// The activity coefficients (natural log) of the species in the system.
    ArrayXr ln_g;

    /// The activities (natural log) of the species in the system.
    ArrayXr ln_a;

    /// The chemical potentials of the species in the system.
    ArrayXr u;

    /// The extra data produced during the evaluation of activity models. This
    /// extra data allows the activity model of a phase to reuse calculated
    /// data from the activity model of a previous phase if needed.
    Map<String, Any> m_extra;
};

/// Output a ChemicalProps object to an output stream.
auto operator<<(std::ostream& out, const ChemicalProps& props) -> std::ostream&;

} // namespace Reaktoro
