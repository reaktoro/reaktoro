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

// C++ includes
#include <string>
#include <map>

// Reaktor includes
#include <Reaktor/Math/BilinearInterpolator.hpp>

namespace Reaktor {

/// A type to describe the common attributes of all kinds of species.
/// The GeneralSpecies is used as a base class for all other species classes.
/// @see AqueousSpecies, GaseousSpecies, MineralSpecies
/// @ingroup Species
struct GeneralSpecies
{
    /// The name of the species
    std::string name;

    /// The chemical formula of the species
    std::string formula;

    /// The elements that compose the species and their number of atoms
    std::map<std::string, double> elements;

    /// The molar mass of the species (in units of kg/mol)
    double molar_mass;

    /// The electrical charge of the species
    double charge;
};

/// A type for storing thermodynamic properties of a species over a range of temperatures and pressures.
/// The thermodynamic properties are represented by BilinearInterpolator instances, where
/// the temperature points (in units of K) are along the x-coordinates, and the pressure points
/// (in units of Pa) are along the y-coordinates.
struct SpeciesThermoProperties
{
    /// The interpolator of the standard molar Gibbs free energy @f$\Delta G_{f}^{\circ}@f$ of the species (in units of J/mol)
    BilinearInterpolator gibbs_energy;

    /// The interpolator of the standard molar Helmholtz free energy @f$\Delta A_{f}^{\circ}@f$ of the species (in units of J/mol)
    BilinearInterpolator helmholtz_energy;

    /// The interpolator of the standard molar internal energy @f$\Delta U_{f}^{\circ}@f$ of the species (in units of J/mol)
    BilinearInterpolator internal_energy;

    /// The interpolator of the standard molar enthalpy @f$\Delta H_{f}^{\circ}@f$ of the species (in units of J/mol)
    BilinearInterpolator enthalpy;

    /// The interpolator of the standard molar entropy @f$ S^{\circ}@f$ of the species (in units of J/K)
    BilinearInterpolator entropy;

    /// The interpolator of the standard molar volume @f$ V^{\circ}@f$ of the species (in units of m3/mol)
    BilinearInterpolator volume;

    /// The interpolator of the standard molar isobaric heat capacity @f$ C_{p}^{\circ}@f$ of the species (in units of J/(mol*K))
    BilinearInterpolator heat_capacity;
};

/// A type for storing thermodynamic properties of a reaction over a range of temperatures and pressures.
/// The thermodynamic properties are represented by BilinearInterpolator instances, where
/// the temperature points (in units of K) are along the x-coordinates, and the pressure points
/// (in units of Pa) are along the y-coordinates.
struct ReactionThermoProperties
{
    /// The equation of the reaction as pairs `(reactant, stoichiometry)`
    std::map<std::string, double> equation;

    /// The interpolator of the equilibrium constant @f$\ln k_{\mathrm{eq}}@f$ of the reaction
    BilinearInterpolator lnk;

    /// The interpolator of the standard molar Gibbs free energy @f$\Delta G_{r}^{\circ}@f$ of the reaction (in units of J/mol)
    BilinearInterpolator gibbs_energy;

    /// The interpolator of the standard molar Helmholtz free energy @f$\Delta A_{r}^{\circ}@f$ of the reaction (in units of J/mol)
    BilinearInterpolator helmholtz_energy;

    /// The interpolator of the standard molar internal energy @f$\Delta U_{r}^{\circ}@f$ of the reaction (in units of J/mol)
    BilinearInterpolator internal_energy;

    /// The interpolator of the standard molar enthalpy @f$\Delta H_{r}^{\circ}@f$ of the reaction (in units of J/mol)
    BilinearInterpolator enthalpy;

    /// The interpolator of the standard molar entropy @f$\Delta S_{r}^{\circ}@f$ of the reaction (in units of J/K)
    BilinearInterpolator entropy;

    /// The interpolator of the standard molar volume @f$\Delta V_{r}^{\circ}@f$ of the reaction (in units of m3/mol)
    BilinearInterpolator volume;

    /// The interpolator of the standard molar isobaric heat capacity @f$ \Delta C_{p}^{\circ}@f$ of the reaction (in units of J/(mol*K))
    BilinearInterpolator heat_capacity;
};

} // namespace Reaktor
