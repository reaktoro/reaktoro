// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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

#include "EquilibriumRestrictions.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Units.hpp>

namespace Reaktoro {
namespace {

/// Compute the amount of a species given a value in mass or in amount.
auto computeSpeciesAmount(const ChemicalSystem& system, Index ispecies, double value, const String& unit)
{
    if(unit == "mol")
        return value;

    if(units::convertible(unit, "kg"))
    {
        const auto molarmass = system.species(ispecies).molarMass(); // in kg/mol
        value = units::convert(value, unit, "kg"); // from some mass unit to kg
        return value / molarmass; // from kg to mol
    }
    else return units::convert(value, unit, "mol"); // from some amount unit to mol
}

} // namespace

EquilibriumRestrictions::EquilibriumRestrictions(const ChemicalSystem& system)
: msystem(system)
{}

auto EquilibriumRestrictions::cannotReact(Index ispecies) -> void
{
    cannotIncrease(ispecies);
    cannotDecrease(ispecies);
}

auto EquilibriumRestrictions::cannotReact(String species) -> void
{
    const auto ispecies = system().species().indexWithName(species);
    cannotReact(ispecies);
}

auto EquilibriumRestrictions::cannotIncrease(Index ispecies) -> void
{
    const auto numspecies = system().species().size();
    errorif(ispecies >= numspecies, "Given species index `", ispecies, "` is out of bounds (number of species is ", numspecies, ").");
    species_cannot_increase_above.erase(ispecies); // remove possible upper bound restriction set before
    species_cannot_increase.insert(ispecies);
}

auto EquilibriumRestrictions::cannotIncrease(String species) -> void
{
    const auto ispecies = system().species().indexWithName(species);
    cannotIncrease(ispecies);
}

auto EquilibriumRestrictions::cannotIncreaseAbove(Index ispecies, double value, String unit) -> void
{
    const auto numspecies = system().species().size();
    errorif(ispecies >= numspecies, "Given species index `", ispecies, "` is out of bounds (number of species is ", numspecies, ").");
    value = computeSpeciesAmount(msystem, ispecies, value, unit);
    species_cannot_increase.erase(ispecies); // remove possible strict upper bound restriction set before
    species_cannot_increase_above.insert_or_assign(ispecies, value);
}

auto EquilibriumRestrictions::cannotIncreaseAbove(String species, double value, String unit) -> void
{
    const auto ispecies = system().species().indexWithName(species);
    cannotIncreaseAbove(ispecies, value, unit);
}

auto EquilibriumRestrictions::cannotDecrease(Index ispecies) -> void
{
    const auto numspecies = system().species().size();
    errorif(ispecies >= numspecies, "Given species index `", ispecies, "` is out of bounds (number of species is ", numspecies, ").");
    species_cannot_decrease_below.erase(ispecies); // remove possible lower bound restriction set before
    species_cannot_decrease.insert(ispecies);
}

auto EquilibriumRestrictions::cannotDecrease(String species) -> void
{
    const auto ispecies = system().species().indexWithName(species);
    cannotDecrease(ispecies);
}

auto EquilibriumRestrictions::cannotDecreaseBelow(Index ispecies, double value, String unit) -> void
{
    const auto numspecies = system().species().size();
    errorif(ispecies >= numspecies, "Given species index `", ispecies, "` is out of bounds (number of species is ", numspecies, ").");
    value = computeSpeciesAmount(msystem, ispecies, value, unit);
    species_cannot_decrease.erase(ispecies); // remove possible strict lower bound restriction set before
    species_cannot_decrease_below.insert_or_assign(ispecies, value);
}

auto EquilibriumRestrictions::cannotDecreaseBelow(String species, double value, String unit) -> void
{
    const auto ispecies = system().species().indexWithName(species);
    cannotDecreaseBelow(ispecies, value, unit);
}

auto EquilibriumRestrictions::canReactFreely(Index ispecies) -> void
{
    canIncreaseFreely(ispecies);
    canDecreaseFreely(ispecies);
}

auto EquilibriumRestrictions::canReactFreely(String species) -> void
{
    const auto ispecies = system().species().indexWithName(species);
    canReactFreely(ispecies);
}

auto EquilibriumRestrictions::canIncreaseFreely(Index ispecies) -> void
{
    const auto numspecies = system().species().size();
    errorif(ispecies >= numspecies, "Given species index `", ispecies, "` is out of bounds (number of species is ", numspecies, ").");
    species_cannot_increase.erase(ispecies);
    species_cannot_increase_above.erase(ispecies);
}

auto EquilibriumRestrictions::canIncreaseFreely(String species) -> void
{
    const auto ispecies = system().species().indexWithName(species);
    canIncreaseFreely(ispecies);
}

auto EquilibriumRestrictions::canDecreaseFreely(Index ispecies) -> void
{
    const auto numspecies = system().species().size();
    errorif(ispecies >= numspecies, "Given species index `", ispecies, "` is out of bounds (number of species is ", numspecies, ").");
    species_cannot_decrease.erase(ispecies);
    species_cannot_decrease_below.erase(ispecies);
}

auto EquilibriumRestrictions::canDecreaseFreely(String species) -> void
{
    const auto ispecies = system().species().indexWithName(species);
    canDecreaseFreely(ispecies);
}

auto EquilibriumRestrictions::system() const -> const ChemicalSystem&
{
    return msystem;
}

auto EquilibriumRestrictions::speciesCannotIncrease() const -> Set<Index> const&
{
    return species_cannot_increase;
}

auto EquilibriumRestrictions::speciesCannotDecrease() const -> Set<Index> const&
{
    return species_cannot_decrease;
}

auto EquilibriumRestrictions::speciesCannotIncreaseAbove() const -> Map<Index, double> const&
{
    return species_cannot_increase_above;
}

auto EquilibriumRestrictions::speciesCannotDecreaseBelow() const -> Map<Index, double> const&
{
    return species_cannot_decrease_below;
}

} // namespace Reaktoro
