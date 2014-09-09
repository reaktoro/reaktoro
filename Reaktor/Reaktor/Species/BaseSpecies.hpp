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

// Reaktor includes
#include <Reaktor/Core/SpeciesComposition.hpp> // todo move to Common

namespace Reaktor {

/// A type to describe the common attributes of all kinds of species
///
/// The BaseSpecies is used as a base class for all other species classes.
///
/// @see AqueousSpecies, GaseousSpecies, MineralSpecies
/// @ingroup Species
struct BaseSpecies
{
	/// The name of the species
    std::string name;

    /// The chemical formula of the species
    std::string formula;

    /// The elemental composition of the species
    SpeciesComposition composition;

    /// The molar mass of the species
    double molar_mass;

    /// The electrical charge of the species
    double charge;
};

} // namespace Reaktor
