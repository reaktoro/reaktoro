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

// C++ includes
#include <iostream>

namespace Reaktor {

// Reaktor forward declarations
class AqueousSpecies;
class GaseousSpecies;
class MineralSpecies;
struct SpeciesThermo;

struct SpeciesThermo
{
	SpeciesThermo();

	/// The standard molar volume @f$ V^{\circ}@f$ of the species (in units of m3/mol)
	double volume;

	/// The standard molar entropy @f$ S^{\circ}@f$ of the species (in units of J/K)
	double entropy;

	/// The apparent standard molar Helmholtz free energy of formation @f$\Delta A_{f}^{\circ}@f$ of the species (in units of J/mol)
	double helmholtz;

	/// The apparent standard molar internal energy of formation @f$\Delta U_{f}^{\circ}@f$ of the species (in units of J/mol)
	double internal_energy;

	/// The apparent standard molar enthalpy of formation @f$\Delta H_{f}^{\circ}@f$ of the species (in units of J/mol)
	double enthalpy;

	/// The apparent standard molar Gibbs free energy of formation @f$\Delta G_{f}^{\circ}@f$ of the species (in units of J/mol)
	double gibbs;

	/// The standard molar isobaric heat capacity @f$ C_{p}^{\circ}@f$ of the species (in units of J/(mol K))
	double cp;
};

/// Output the thermodynamic state of the species
auto operator<<(std::ostream& out, const SpeciesThermo& st) -> std::ostream&;

/// Calculate the thermodynamic state of the aqueous species using the HKF model
auto speciesThermo(double T, double P, const AqueousSpecies& species) -> SpeciesThermo;

/// Calculate the thermodynamic state of the gaseous species using the HKF model
auto speciesThermo(double T, double P, const GaseousSpecies& species) -> SpeciesThermo;

/// Calculate the thermodynamic state of the mineral species using the HKF model
auto speciesThermo(double T, double P, const MineralSpecies& species) -> SpeciesThermo;

} /* namespace Reaktor */
