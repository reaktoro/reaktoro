/*
 * Reaktor is a C++ library for computational reaction modelling.
 *
 * Copyright (C) 2014 Allan Leal
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

namespace Reaktor {

// Forward declarations
class AqueousSpecies;
struct FunctionG;

struct SpeciesElectroState
{
    SpeciesElectroState();

    /// The effective electrostatic radius of the solute species at referente temperature 298.15 K and pressure 1 bar
    double reref;

    /// The effective electrostatic radius of the solute species
    double re;

    /// The Born coefficient of the solute species
    double w;

    /// The first-order partial derivative of the Born coefficient of the solute species with respect to temperature
    double wT;

    /// The first-order partial derivative of the Born coefficient of the solute species with respect to pressure
    double wP;

    /// The second-order partial derivative of the Born coefficient of the solute species with respect to temperature
    double wTT;

    /// The second-order partial derivative of the Born coefficient of the solute species with respect to temperature and pressure
    double wTP;

    /// The second-order partial derivative of the Born coefficient of the solute species with respect to pressure
    double wPP;
};

/// Calculate the electrostatic state of the aqueous species
auto speciesElectro(double T, double P, const AqueousSpecies& species) -> SpeciesElectroState;

/// Calculate the electrostatic state of the aqueous species using the g-function state
auto speciesElectro(const FunctionG& g, const AqueousSpecies& species) -> SpeciesElectroState;

} // namespace Reaktor
