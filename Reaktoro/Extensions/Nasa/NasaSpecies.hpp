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
#include <Reaktoro/Extensions/Nasa/NasaThermoData.hpp>

namespace Reaktoro {

// Forward declarations
class Species;

/// The possible aggregate states of a species in a NASA thermodynamic database.
enum class NasaAggregateState { Gas, Liquid, Solid };

/// The possible types of species in a NASA thermodynamic database.
enum class NasaSpeciesType { Product, Reactant, ProductReactant };

/// Used to represent a species in a NASA thermodynamic database.
struct NasaSpecies
{
    /// The name of the species.
    String name;

    /// The comment and data source of the species.
    String comment;

    /// The identification code of the species.
    String idcode;

    /// The element symbols and coefficients in the formula of the species.
    Pairs<String, double> formula;

    /// The aggregate state code of the species (zero for gas, nonzero for condensed phases).
    int aggregatecode;

    /// The aggregate state of the species (liquid, gas, or solid) determined from name and aggregate state code.
    NasaAggregateState aggregatestate;

    /// The type of the species.
    NasaSpeciesType type;

    /// The molar mass of the species (in g/mol).
    real molarmass;

    /// The heat of formation \eq{\Delta H_{f}^{\circ}} at 298.15 K (in J/mol).
    real dHf;

    /// The value of \eq{\Delta H_{0}^{\circ}=H^{\circ}(298.15)-H^{\circ}(0)} (in J/mol).
    real dH0;

    /// The assigned enthalpy (in J/mol) of the species when there are no temperature intervals.
    real H0;

    /// The temperature (in K) corresponding to the assigned enthalpy when there are no temperature intervals.
    real TH0;

    /// The data used to compute standard thermodynamic properties of the species at different temperature ranges.
    Vec<NasaThermoData> thermodata;
};

/// Used to represent a species in a NASA thermodynamic database that has or not phase transitions.
struct NasaExtendedSpecies
{
    /// The original name of the species if it has no phase transition, otherwise a name with suffix (cd) to denote condensed substance.
    String name;

    /// The different species to represent this extended species over a wider range of temperature.
    Vec<NasaSpecies> specieslist;

    /// Return the aggregate state of the extended species at given temperature.
    auto aggregateStateAtTemperature(double T) const -> NasaAggregateState;

    /// Return the aggregate state of the extended species at given temperature.
    auto speciesNameAtTemperature(double T) const -> String;
};

/// Return true if two NasaSpecies objects are different.
auto operator!=(const NasaSpecies& l, const NasaSpecies& r) -> bool;

/// Return true if two NasaSpecies objects are equal.
auto operator==(const NasaSpecies& l, const NasaSpecies& r) -> bool;

} // namespace Reaktoro
