// Reaktoro is a unified framework for modeling chemically reactive phases.
//
// Copyright © 2014-2024 Allan Leal
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
#include <Reaktoro/Common/Types.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalProps;
class ChemicalSystem;

/// The function type for evaluation of a property of a chemical system.
/// @param props The already evaluated chemical properties of the system.
/// @return The property of interest retrieved from the given chemical properties of the system.
using PropFn = Fn<real(const ChemicalProps& props)>;

/// Used to retrieve a specified property from the chemical properties of a system.
class Prop
{
public:
    /// Construct a Prop evaluator object with given property evaluation function.
    explicit Prop(const PropFn& propfn);

    /// Evaluate the property with given chemical properties of the system.
    /// @param props The already evaluated chemical properties of the system.
    auto eval(const ChemicalProps& props) -> real;

    /// Evaluate the property with given chemical properties of the system.
    /// @param props The already evaluated chemical properties of the system.
    auto operator()(const ChemicalProps& props) -> real;

    /// Return a property function that evaluates the amount of an element in the system (in mol).
    static auto elementAmount(const ChemicalSystem& system, const String& element) -> Prop;

    /// Return a property function that evaluates the amount of an element in a phase of the system (in mol).
    static auto elementAmountInPhase(const ChemicalSystem& system, const String& element, const String& phase) -> Prop;

    /// Return a property function that evaluates the mass of an element in the system (in kg).
    static auto elementMass(const ChemicalSystem& system, const String& element) -> Prop;

    /// Return a property function that evaluates the mass of an element in a phase of the system (in kg).
    static auto elementMassInPhase(const ChemicalSystem& system, const String& element, const String& phase) -> Prop;

    /// Return a property function that evaluates the amount of a species in the system (in mol).
    static auto speciesAmount(const ChemicalSystem& system, const String& species) -> Prop;

    /// Return a property function that evaluates the mass of a species in the system (in kg).
    static auto speciesMass(const ChemicalSystem& system, const String& species) -> Prop;

    /// Return a property function that evaluates the mole fraction of a species in the system.
    static auto speciesMoleFraction(const ChemicalSystem& system, const String& species) -> Prop;

    /// Return a property function that evaluates the activity coefficient of a species in the system.
    static auto speciesActivityCoefficient(const ChemicalSystem& system, const String& species) -> Prop;

    /// Return a property function that evaluates the activity of a species in the system.
    static auto speciesActivity(const ChemicalSystem& system, const String& species) -> Prop;

    /// Return a property function that evaluates the chemical potentials of a species in the system (in J/mol).
    static auto speciesChemicalPotential(const ChemicalSystem& system, const String& species) -> Prop;

    /// Return a property function that evaluates the standard partial molar volume of a species in the system (in m3/mol).
    static auto speciesStandardVolume(const ChemicalSystem& system, const String& species) -> Prop;

    /// Return a property function that evaluates the standard partial molar Gibbs energy of formation of a species in the system (in J/mol).
    static auto speciesStandardGibbsEnergy(const ChemicalSystem& system, const String& species) -> Prop;

    /// Return a property function that evaluates the standard partial molar enthalpy of formation of a species in the system (in J/mol).
    static auto speciesStandardEnthalpy(const ChemicalSystem& system, const String& species) -> Prop;

    /// Return a property function that evaluates the standard partial molar entropy of formation of a species in the system (in J/(mol*K)).
    static auto speciesStandardEntropy(const ChemicalSystem& system, const String& species) -> Prop;

    /// Return a property function that evaluates the standard partial molar internal energy of formation of a species in the system (in J/mol).
    static auto speciesStandardInternalEnergy(const ChemicalSystem& system, const String& species) -> Prop;

    /// Return a property function that evaluates the standard partial molar Helmholtz energy of formation of a species in the system (in J/mol).
    static auto speciesStandardHelmholtzEnergy(const ChemicalSystem& system, const String& species) -> Prop;

    /// Return a property function that evaluates the standard partial molar isobaric heat capacity of a species in the system (in J/(mol*K)).
    static auto speciesStandardHeatCapacitiesConstP(const ChemicalSystem& system, const String& species) -> Prop;

    /// Return a property function that evaluates the standard partial molar isochoric heat capacity of a species in the system (in J/(mol*K)).
    static auto speciesStandardHeatCapacitiesConstV(const ChemicalSystem& system, const String& species) -> Prop;

    /// Return a property function that evaluates the sum of species amounts in a phase of the system (in mol).
    static auto phaseAmount(const ChemicalSystem& system, const String& phase) -> Prop;

    /// Return a property function that evaluates the sum of species masses in a phase of the system (in kg).
    static auto phaseMass(const ChemicalSystem& system, const String& phase) -> Prop;

    /// Return a property function that evaluates the volume of a phase in the system (in m3).
    static auto phaseVolume(const ChemicalSystem& system, const String& phase) -> Prop;

    /// Return a property function that evaluates the Gibbs energy of formation of a phase in the system (in J).
    static auto phaseGibbsEnergy(const ChemicalSystem& system, const String& phase) -> Prop;

    /// Return a property function that evaluates the enthalpy of formation of a phase in the system (in J).
    static auto phaseEnthalpy(const ChemicalSystem& system, const String& phase) -> Prop;

    /// Return a property function that evaluates the entropy of formation of the a phase in system (in J/K).
    static auto phaseEntropy(const ChemicalSystem& system, const String& phase) -> Prop;

    /// Return a property function that evaluates the internal energy of formation of a phase in the system (in J).
    static auto phaseInternalEnergy(const ChemicalSystem& system, const String& phase) -> Prop;

    /// Return a property function that evaluates the Helmholtz energy of formation of a phase in the system (in J).
    static auto phaseHelmholtzEnergy(const ChemicalSystem& system, const String& phase) -> Prop;

    /// Return a property function that evaluates the temperature of the system (in K).
    static auto temperature(const ChemicalSystem& system) -> Prop;

    /// Return a property function that evaluates the pressure of the system (in Pa).
    static auto pressure(const ChemicalSystem& system) -> Prop;

    /// Return a property function that evaluates the sum of species amounts in the system (in mol).
    static auto amount(const ChemicalSystem& system) -> Prop;

    /// Return a property function that evaluates the sum of species masses in the system (in kg).
    static auto mass(const ChemicalSystem& system) -> Prop;

    /// Return a property function that evaluates the volume of the system (in m3).
    static auto volume(const ChemicalSystem& system) -> Prop;

    /// Return a property function that evaluates the Gibbs energy of formation of the system (in J).
    static auto gibbsEnergy(const ChemicalSystem& system) -> Prop;

    /// Return a property function that evaluates the enthalpy of formation of the system (in J).
    static auto enthalpy(const ChemicalSystem& system) -> Prop;

    /// Return a property function that evaluates the entropy of formation of the system (in J/K).
    static auto entropy(const ChemicalSystem& system) -> Prop;

    /// Return a property function that evaluates the internal energy of formation of the system (in J).
    static auto internalEnergy(const ChemicalSystem& system) -> Prop;

    /// Return a property function that evaluates the Helmholtz energy of formation of the system (in J).
    static auto helmholtzEnergy(const ChemicalSystem& system) -> Prop;

private:
    /// The function that evaluates/retrieves a property of a chemical system.
    PropFn propfn;
};

} // namespace Reaktoro
