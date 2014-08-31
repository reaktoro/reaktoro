// Reaktor is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
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

// Reaktor includes
#include <Reaktor/Thermodynamics/ThermoDataReaction.hpp>

namespace Reaktor {

// Forward declarations
class BilinearInterpolator;

/// A type for storing thermodynamic properties of a species over a range of temperatures and pressures
///
/// The thermodynamic properties are represented by BilinearInterpolator instances, where
/// the temperature points (in units of K) are along the x-coordinates, and the pressure points
/// (in units of Pa) are along the y-coordinates.
///
/// @see BilinearInterpolator, ThermoDataReaction
struct ThermoDataSpecies
{
    /// The interpolator of the apparent standard molar Gibbs free energy of formation
    /// @f$\Delta G_{f}^{\circ}@f$ of the species (in units of J/mol)
    BilinearInterpolator gibbs;

    /// The interpolator of the apparent standard molar enthalpy of formation
    /// @f$\Delta H_{f}^{\circ}@f$ of the species (in units of J/mol)
    BilinearInterpolator enthalpy;

    /// The interpolator of the standard molar entropy
    /// @f$ S^{\circ}@f$ of the species (in units of J/K)
    BilinearInterpolator entropy;

    /// The interpolator of the standard molar volume
    /// @f$ V^{\circ}@f$ of the species (in units of m3/mol)
    BilinearInterpolator volume;

    /// The interpolator of the standard molar isobaric heat capacity
    /// @f$ C_{p}^{\circ}@f$ of the species (in units of J/(mol*K))
    BilinearInterpolator cp;

    /// The thermodynamic properties of the species given in terms of a reaction
    ThermoDataReaction reaction;
};

} // namespace Reaktor
