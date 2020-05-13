// Reaktoro is a unified framework for modeling chemically reactive phases.
//
// Copyright (C) 2014-2020 Allan Leal
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
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/PhaseThermoProps.hpp>

namespace Reaktoro {

/// The primary standard thermodynamic property data of the phases and species in a chemical system.
/// @see ThermoProps
struct ThermoPropsData
{
    /// The temperature of the system (in K).
    real T = {};

    /// The pressure of the system (in Pa).
    real P = {};

    /// The standard molar Gibbs energies of the species in the system (in J/mol).
    ArrayXr G0;

    /// The standard molar enthalpies of the species in the system (in J/mol).
    ArrayXr H0;

    /// The standard molar volumes of the species in the system (in m3/mol).
    ArrayXr V0;

    /// The standard molar isobaric heat capacities of the species in the system (in J/(mol·K)).
    ArrayXr Cp0;

    /// The standard molar isochoric heat capacities of the species in the system (in J/(mol·K)).
    ArrayXr Cv0;
};

/// The standard thermodynamic properties of the phases and species in a chemical system.
class ThermoProps
{
public:
    /// Construct a ThermoProps object.
    explicit ThermoProps(const ChemicalSystem& system);

    /// Construct a ThermoProps object.
    ThermoProps(const ChemicalSystem& system, const ThermoPropsData& data);

    /// Update the standard thermodynamic properties of the phases and its species in the chemical system.
    /// @param T The temperature condition (in K)
    /// @param P The pressure condition (in Pa)
    auto update(real T, real P) -> void;

    /// Return the chemical system associated with these standard thermodynamic properties.
    auto system() const -> const ChemicalSystem&;

    /// Return the primary standard thermodynamic property data of the system from which others are calculated.
    auto data() const -> const ThermoPropsData&;

    /// Return the standard thermodynamic properties of a phase with given index.
    auto phaseProps(Index idx) const -> PhaseThermoPropsConstRef;

    /// Return the standard thermodynamic properties of a phase with given index.
    auto phaseProps(Index idx) -> PhaseThermoPropsRef;

    /// Return the temperature of the system (in K).
    auto temperature() const -> real;

    /// Return the pressure of the system (in Pa).
    auto pressure() const -> real;

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

private:
    /// The chemical system associated with these standard thermodynamic properties.
    ChemicalSystem sys;

    /// The primary standard thermodynamic property data of the system from which others are calculated.
    ThermoPropsData props;
};

} // namespace Reaktoro
