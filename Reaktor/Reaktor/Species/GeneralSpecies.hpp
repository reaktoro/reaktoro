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
#include <vector>

namespace Reaktor {

/// Define a base class for all species classes
///
/// The GeneralSpecies class is used to represent the common properties of
/// chemical species. It is used as a base class for all other species classes.
///
/// @see AqueousSpecies, GaseousSpecies, MineralSpecies
/// @ingroup Species
struct GeneralSpecies
{
	/// Construct a default GeneralSpecies instance
	GeneralSpecies();

	/// The name of the chemical species
    std::string name;

    /// The formula of the chemical species
    std::string formula;

    /// The elements that compose the chemical species
    std::vector<std::string> elements;

    /// The coefficients of the elements that compose the chemical species
    std::vector<double> coefficients;

    /// The molar mass of the chemical species
    double molar_mass;

    /// The electrical charge of the chemical species
    double charge;
};

} // namespace Reaktor
