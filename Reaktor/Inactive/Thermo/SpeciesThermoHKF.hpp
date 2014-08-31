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
class GaseousSpecies;
class MineralSpecies;
struct SpeciesElectro;
struct SpeciesThermo;
struct WaterElectro;
struct WaterThermo;

/// Calculate the thermodynamic state of the aqueous species, but not H2O(l), using the HKF model
auto speciesThermoHKF(double T, double P, const AqueousSpecies& species, const SpeciesElectro& se, const WaterElectro& we) -> SpeciesThermo;

/// Calculate the thermodynamic state of the aqueous species using the HKF model
auto speciesThermoHKF(double T, double P, const AqueousSpecies& species) -> SpeciesThermo;

/// Calculate the thermodynamic state of the gaseous species using the HKF model
auto speciesThermoHKF(double T, double P, const GaseousSpecies& species) -> SpeciesThermo;

/// Calculate the thermodynamic state of the mineral species using the HKF model
auto speciesThermoHKF(double T, double P, const MineralSpecies& species) -> SpeciesThermo;

/// Calculate the thermodynamic state of the water species H2O(l) using the HKF model and the Wagner and Pruss (1995) equation of state
auto speciesThermoHKF(double T, double P, const WaterThermo& wt) -> SpeciesThermo;

} // namespace Reaktor
