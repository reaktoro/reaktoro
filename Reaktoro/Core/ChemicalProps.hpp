// Reaktoro is a unified framework for modeling chemically reactive phases.
//
// Copyright (C) 2014-2021 Allan Leal
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
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ChemicalPropsPhase.hpp>

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

    /// Construct a copy of a ChemicalProps object.
    ChemicalProps(const ChemicalProps& other);

    /// Destroy this ChemicalProps object.
    ~ChemicalProps();

    /// Assign a ChemicalProps object to this.
    auto operator=(ChemicalProps other) -> ChemicalProps&;

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

    /// Update the chemical properties of the system using the array stream @p stream.
    /// @param stream The array stream containing the serialized chemical properties.
    auto deserialize(const ArrayStream<real>& stream) -> void;

    /// Return the chemical system associated with these chemical properties.
    auto system() const -> const ChemicalSystem&;

    /// Return the chemical properties of a phase with given index.
    auto phaseProps(Index idx) const -> ChemicalPropsPhaseConstRef;

    /// Return the temperature of the system (in K).
    auto temperature() const -> const real&;

    /// Return the pressure of the system (in Pa).
    auto pressure() const -> const real&;

    /// Return the amounts of each species in the system (in mol).
    auto speciesAmounts() const -> ArrayXrConstRef;

    /// Return the mole fractions of the species in the system.
    auto moleFractions() const -> ArrayXrConstRef;

    /// Return the ln activity coefficients of the species in the system.
    auto lnActivityCoefficients() const -> ArrayXrConstRef;

    /// Return the ln activities of the species in the system.
    auto lnActivities() const -> ArrayXrConstRef;

    /// Return the chemical potentials of the species in the system (in J/mol).
    auto chemicalPotentials() const -> ArrayXrConstRef;

    /// Return the standard partial molar Gibbs energies of the species in the system (in J/mol).
    auto standardGibbsEnergies() const -> ArrayXrConstRef;

    /// Return the standard partial molar enthalpies of the species in the system (in J/mol).
    auto standardEnthalpies() const -> ArrayXrConstRef;

    /// Return the standard partial molar volumes of the species in the system (in m3/mol).
    auto standardVolumes() const -> ArrayXrConstRef;

    /// Return the standard partial molar entropies of the species in the system (in J/(mol*K)).
    auto standardEntropies() const -> ArrayXr;

    /// Return the standard partial molar internal energies of the species in the system (in J/mol).
    auto standardInternalEnergies() const -> ArrayXr;

    /// Return the standard partial molar Helmholtz energies of the species in the system (in J/mol).
    auto standardHelmholtzEnergies() const -> ArrayXr;

    /// Return the standard partial molar isobaric heat capacities of the species in the system (in J/(mol*K)).
    auto standardHeatCapacitiesConstP() const -> ArrayXrConstRef;

    /// Return the standard partial molar isochoric heat capacities of the species in the system (in J/(mol*K)).
    auto standardHeatCapacitiesConstV() const -> ArrayXrConstRef;

    /// Return the Gibbs energy of the system (in J).
    auto gibbsEnergy() const -> real;

    /// Return the enthalpy of the system (in J).
    auto enthalpy() const -> real;

    /// Return the volume of the system (in m3).
    auto volume() const -> real;

    /// Return the entropy of the system (in J/K).
    auto entropy() const -> real;

    /// Return the internal energy of the system (in J).
    auto internalEnergy() const -> real;

    /// Return the Helmholtz energy of the system (in J).
    auto helmholtzEnergy() const -> real;

    /// Return the chemical properties in this object serialized in an array of real numbers.
    operator VectorXr() const;

    /// Return the chemical properties in this object serialized in an array of double numbers.
    operator VectorXd() const;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

} // namespace Reaktoro
