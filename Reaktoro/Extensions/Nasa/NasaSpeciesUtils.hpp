// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/AggregateState.hpp>

namespace Reaktoro {

// Forward declarations
class Element;
class ElementalComposition;
class NasaSpecies;
class Species;

namespace NasaUtils {

/// Correct all uppercase element symbols in NASA database to their correct name.
auto correctElementSymbol(const String& symbol) -> String;

/// Return an Element object for given element symbol used in the NASA database.
auto createElement(const String& symbol) -> Element;

/// Return the elemental composition of a chemical species from the NASA database.
auto createElements(const NasaSpecies& species) -> ElementalComposition;

/// Return the electric charge of a chemical species from the NASA database.
auto charge(const NasaSpecies& species) -> double;

/// Return the aggregate state of a chemical species from the NASA database.
auto aggregateState(const NasaSpecies& species) -> AggregateState;

/// Return the tags for a chemical species from the NASA database.
/// This method returns tags `{"product"}` for a product species and
/// `{"reactant"}` for a reactant species.
auto tags(const NasaSpecies& species) -> Strings;

/// Convert this NasaSpecies object into a Species object.
auto convertSpecies(const NasaSpecies& species) -> Species;

} // namespace NasaUtils
} // namespace Reaktoro
