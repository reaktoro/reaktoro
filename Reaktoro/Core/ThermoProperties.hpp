// Reaktoro is a unified framework for modeling chemically reactive systems.
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

// C++ includes
#include <memory>

// Reaktoro includes
#include <Reaktoro/Common/ThermoVector.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalSystem;

/// A class used for calculating standard thermodynamic properties of species in a chemical system.
class ThermoProperties
{
public:
    /// Construct a default ThermoProperties instance.
    ThermoProperties();

    /// Construct a ThermoProperties instance with given ChemicalSystem.
    ThermoProperties(const ChemicalSystem& system);

    /// Construct a copy of a ThermoProperties instance.
    ThermoProperties(const ThermoProperties& other);

    /// Destroy this instance.
    virtual ~ThermoProperties();

    /// Construct a copy of a ThermoProperties instance.
    auto operator=(ThermoProperties other) -> ThermoProperties&;

    /// Update the thermodynamic properties of the chemical system.
    /// @param T The new temperature (in units of K)
    /// @param P The new pressure (in units of Pa)
    auto update(double T, double P) -> void;

    /// Return the temperature of the phase (in units of K).
    auto temperature() const -> double;

    /// Return the pressure of the phase (in units of Pa).
    auto pressure() const -> double;

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

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
