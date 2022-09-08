// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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
#include <Reaktoro/Core/ReactionStandardThermoProps.hpp>
#include <Reaktoro/Extensions/Phreeqc/PhreeqcLegacy.hpp> // ***BECAUSE PhreeqcLegacy.hpp IS INCLUDED HERE, MAKE SURE THIS HEADER FILE IS NOT EXPORTED!***
#include <Reaktoro/Extensions/Phreeqc/PhreeqcWater.hpp>

namespace Reaktoro {
namespace PhreeqcUtils {

/// Return the standard molar volume of a PHREEQC species using same model as used in PHREEQC (in cm3/mol).
/// @param species The pointer to the Phreeqc species
/// @param T The temperature value (in K)
/// @param P The pressure value (in Pa)
/// @param wprops The thermodynamic and electrostatic properties of water
auto standardVolume(const PhreeqcSpecies* species, real T, real P, const PhreeqcWaterProps& wprops) -> real;

/// Return the standard molar volume of a PHREEQC phase using same model as used in PHREEQC (in cm3/mol).
/// @param phase The pointer to the Phreeqc phase
/// @param T The temperature value (in K)
/// @param P The pressure value (in Pa)
auto standardVolume(const PhreeqcPhase* phase, real T, real P) -> real;

/// Create a standard thermodynamic model for a PHREEQC species.
auto reactionThermoModel(const PhreeqcSpecies* species) -> ReactionStandardThermoModel;

/// Create a standard thermodynamic model for a PHREEQC phase.
auto reactionThermoModel(const PhreeqcPhase* phase) -> ReactionStandardThermoModel;

/// Create a standard molar volume model for a PHREEQC species.
auto standardVolumeModel(const PhreeqcSpecies* species) -> Model<real(real,real)>;

/// Create a standard molar volume model for a PHREEQC phase.
auto standardVolumeModel(const PhreeqcPhase* phase) -> Model<real(real,real)>;

} // namespace PhreeqcUtils
} // namespace Reaktoro
