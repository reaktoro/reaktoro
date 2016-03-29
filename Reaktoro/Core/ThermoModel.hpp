// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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
#include <Reaktoro/Common/ThermoVector.hpp>

namespace Reaktoro {

// Forward declarations
class PhaseThermoModelResult;

/// The result of the thermodynamic model function of a chemical system.
/// This class holds the standard thermodynamic properties of the species
/// in a chemical system calculated by a ThermoModel function.
/// @see ThermoModel, ChemicalSystem
struct ThermoModelResult
{
    /// Construct a default ThermoModelResult instance
    ThermoModelResult();

    /// Construct a ThermoModelResult instance with allocated memory
    explicit ThermoModelResult(unsigned nspecies);

    /// Assign a vector of PhaseThermoModelResult instances to this.
    auto operator=(const std::vector<PhaseThermoModelResult>& results) -> ThermoModelResult&;

    /// Return the standard partial molar Gibbs energies of the species (in units of J/mol).
    auto standardPartialMolarGibbsEnergies() const -> ThermoVector;

    /// Return the standard partial molar enthalpies of the species (in units of J/mol).
    auto standardPartialMolarEnthalpies() const -> ThermoVector;

    /// Return the standard partial molar volumes of the species (in units of m3/mol).
    auto standardPartialMolarVolumes() const -> ThermoVector;

    /// Return the standard partial molar isobaric heat capacities of the species (in units of J/(mol*K)).
    auto standardPartialMolarHeatCapacitiesConstP() const -> ThermoVector;

    /// Return the standard partial molar isochoric heat capacities of the species (in units of J/(mol*K)).
    auto standardPartialMolarHeatCapacitiesConstV() const -> ThermoVector;

private:
    /// The number of species in the system.
    unsigned num_species = 0;

    /// The results of the phase thermo model function of each phase.
    std::vector<PhaseThermoModelResult> results;
};

} // namespace Reaktoro
