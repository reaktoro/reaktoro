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

namespace Reaktor {

// Forward declarations
class BilinearInterpolator;

/// A type to describe the equation of a reaction // todo move this to somewhere else
struct ReactionEquation
{
	/// The names of the reactants in the reaction
    std::vector<std::string> reactants;

    /// The stoichiometries of the reactants in the reaction
    std::vector<double> stoichiometries;
};

/// A type for storing thermodynamic properties of a reaction over a range of temperatures and pressures
///
/// The thermodynamic properties are represented by BilinearInterpolator instances, where
/// the temperature points (in units of K) are along the x-coordinates, and the pressure points
/// (in units of Pa) are along the y-coordinates.
///
/// @see BilinearInterpolator, ThermoDataSpecies
struct ThermoDataReaction
{
	/// The equation of the reaction with its reactants and stoichiometries
    ReactionEquation equation;

    /// The interpolator of the equilibrium constant @f$\log K_{\mathrm{eq}}@f$ of the reaction
    BilinearInterpolator logK;

    /// The interpolator of the equilibrium constant @f$\ln K_{\mathrm{eq}}@f$ of the reaction
    BilinearInterpolator lnK;

    /// The interpolator of the equilibrium constant @f$\mathrm{p}K_{\mathrm{eq}}@f$ of the reaction
    BilinearInterpolator pK;

    /// The interpolator of the standard molar Gibbs free energy @f$\Delta G_{r}^{\circ}@f$ of the reaction (in units of J/mol)
    BilinearInterpolator gibbs;

    /// The interpolator of the standard molar enthalpy @f$\Delta H_{r}^{\circ}@f$ of the reaction (in units of J/mol)
    BilinearInterpolator enthalpy;

    /// The interpolator of the standard molar entropy @f$\Delta S_{r}^{\circ}@f$ of the reaction (in units of J/K)
    BilinearInterpolator entropy;

    /// The interpolator of the standard molar volume @f$\Delta V_{r}^{\circ}@f$ of the reaction (in units of m3/mol)
    BilinearInterpolator volume;
};

} // namespace Reaktor
