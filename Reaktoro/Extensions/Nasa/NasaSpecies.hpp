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
#include <Reaktoro/Extensions/Nasa/NasaThermoParams.hpp>

namespace Reaktoro {

// Forward declarations
class Species;

/// The possible aggregate states of a chemical species in a NASA thermodynamic database.
enum class NasaAggregateState { Gas, Condensed };

/// The possible types of chemical species in a NASA thermodynamic database.
enum class NasaSpeciesType { Product, Reactant, ProductReactant };

/// Used to represent a species in a NASA thermodynamic database.
struct NasaSpecies
{
    /// The name of the chemical species.
    String name;

    /// The comment and data source of the chemical species.
    String comment;

    /// The identification code of the chemical species.
    String idcode;

    /// The element symbols and coefficients in the formula of the chemical species.
    Pairs<String, double> formula;

    /// The aggregate state of the chemical species.
    NasaAggregateState aggregatestate;

    /// The type of the chemical species.
    NasaSpeciesType type;

    /// The molar mass of the species (in g/mol).
    real molarmass;

    /// The heat of formation \eq{\Delta H_{f}^{\circ}} at 298.15 K (in J/mol).
    real dHf;

    /// The value of \eq{\Delta H_{0}^{\circ}=H^{\circ}(298.15)-H^{\circ}(0)} (in J/mol).
    real dH0;

    /// The assigned enthalpy (in J/mol) of the chemical species when there are no temperature intervals.
    real H0;

    /// The temperature (in K) corresponding to the assigned enthalpy when there are no temperature intervals.
    real TH0;

    /// The data used to compute standard thermodynamic properties of the species at different temperature ranges.
    Vec<NasaThermoParams> thermodata;
};

/// Return true if two NasaSpecies objects are different.
auto operator!=(const NasaSpecies& l, const NasaSpecies& r) -> bool;

/// Return true if two NasaSpecies objects are equal.
auto operator==(const NasaSpecies& l, const NasaSpecies& r) -> bool;

} // namespace Reaktoro
