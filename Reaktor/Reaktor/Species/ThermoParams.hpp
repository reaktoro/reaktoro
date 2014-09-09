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
#include <Reaktor/Common/Optional.hpp>
#include <Reaktor/Core/ReactionEquation.hpp>
#include <Reaktor/Math/BilinearInterpolator.hpp>

namespace Reaktor {

/// A type for storing thermodynamic properties of a reaction over a range of temperatures and pressures
///
/// The thermodynamic properties are represented by BilinearInterpolator instances, where
/// the temperature points (in units of K) are along the x-coordinates, and the pressure points
/// (in units of Pa) are along the y-coordinates.
struct ThermoParamsReaction
{
    /// The equation of the reaction
    Optional<ReactionEquation> equation;

    /// The interpolator of the equilibrium constant @f$\log K_{\mathrm{eq}}@f$ of the reaction
    Optional<BilinearInterpolator> logK;

    /// The interpolator of the equilibrium constant @f$\ln K_{\mathrm{eq}}@f$ of the reaction
    Optional<BilinearInterpolator> lnK;

    /// The interpolator of the equilibrium constant @f$\mathrm{p}K_{\mathrm{eq}}@f$ of the reaction
    Optional<BilinearInterpolator> pK;

    /// The interpolator of the standard molar Gibbs free energy @f$\Delta G_{r}^{\circ}@f$ of the reaction (in units of J/mol)
    Optional<BilinearInterpolator> gibbs;

    /// The interpolator of the standard molar enthalpy @f$\Delta H_{r}^{\circ}@f$ of the reaction (in units of J/mol)
    Optional<BilinearInterpolator> enthalpy;

    /// The interpolator of the standard molar entropy @f$\Delta S_{r}^{\circ}@f$ of the reaction (in units of J/K)
    Optional<BilinearInterpolator> entropy;

    /// The interpolator of the standard molar volume @f$\Delta V_{r}^{\circ}@f$ of the reaction (in units of m3/mol)
    Optional<BilinearInterpolator> volume;
};

/// A type for storing thermodynamic properties of a species over a range of temperatures and pressures
///
/// The thermodynamic properties are represented by BilinearInterpolator instances, where
/// the temperature points (in units of K) are along the x-coordinates, and the pressure points
/// (in units of Pa) are along the y-coordinates.
struct ThermoParamsSpecies
{
	/// The interpolator of the standard molar Gibbs energy @f$\Delta G_{f}^{\circ}@f$ of the species (in units of J/mol)
	Optional<BilinearInterpolator> gibbs;

	/// The interpolator of the standard molar enthalpy @f$\Delta H_{f}^{\circ}@f$ of the species (in units of J/mol)
	Optional<BilinearInterpolator> enthalpy;

	/// The interpolator of the standard molar entropy @f$ S^{\circ}@f$ of the species (in units of J/K)
	Optional<BilinearInterpolator> entropy;

	/// The interpolator of the standard molar volume @f$ V^{\circ}@f$ of the species (in units of m3/mol)
	Optional<BilinearInterpolator> volume;

	/// The interpolator of the standard molar isobaric heat capacity @f$ C_{p}^{\circ}@f$ of the species (in units of J/(mol*K))
	Optional<BilinearInterpolator> cp;

	/// The thermodynamic properties of the species in terms of the thermodynamic properties of a reaction
	ThermoParamsReaction reaction;
};

} // namespace Reaktor
